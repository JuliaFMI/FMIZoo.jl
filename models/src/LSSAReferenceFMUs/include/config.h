#ifndef config_h
#define config_h

// define class name and unique id
// #define MODEL_IDENTIFIER BouncingBall
#define INSTANTIATION_TOKEN "{1AE5E10D-9521-4DE3-80B9-D0EAAA7D5AF1}"

#define CO_SIMULATION
#define MODEL_EXCHANGE

#define HAS_CONTINUOUS_STATES
#define HAS_EVENT_INDICATORS

#define SET_FLOAT64
#define GET_OUTPUT_DERIVATIVE
#define EVENT_UPDATE

#define FIXED_SOLVER_STEP 1e-3
#define DEFAULT_STOP_TIME 3

#ifndef GET_PARTIAL_DERIVATIVE
#define GET_PARTIAL_DERIVATIVE
#endif

// LS-SA: sparse fmi3GetDirectionalDerivative / fmi3GetAdjointDerivative implemented in model.c
#define LSSA_DERIVATIVES

#define GET_INT32
#define SET_INT32

#ifndef NUM_BALLS
#define NUM_BALLS 1
#endif

// Global Variables
typedef enum {
    vr_time = 0,
    vr_g = 1,
    vr_v_min = 2,
    vr_ground_flexibility = 3,
    vr_base_event_period = 4,
    vr_nextEventTime = 5
} ValueReferenceGlobal;

// number of global value references (0 ... NUM_GLOBAL_VRS-1)
#define NUM_GLOBAL_VRS 6

#ifndef START_VR
#define START_VR 100
#endif

// Coupled variant ("NCoupledBouncingBalls", NUM_NEIGHBORS > 0): every ball is connected to its NUM_NEIGHBORS
// nearest neighbors (by index) on each side by a hardening spring with damper, and passes a share of its
// impact impulse to them. NUM_NEIGHBORS = 0 is the uncoupled model (identical code and variables).
#ifndef NUM_NEIGHBORS
#define NUM_NEIGHBORS 0
#endif

#if NUM_NEIGHBORS > 0
#define VARS_PER_BALL 13
#else
#define VARS_PER_BALL 12
#endif

#define OFF_H 0
#define OFF_DER_H 1
#define OFF_V 2
#define OFF_DER_V 3
#define OFF_E 4
#define OFF_K 5
#define OFF_DRAG 6
#define OFF_BOUNCE 7
#define OFF_PRE_H 8
#define OFF_PRE_V 9
#define OFF_ERR_H 10
#define OFF_EVENT_G 11
#if NUM_NEIGHBORS > 0
#define OFF_C 12            // spring stiffness c_i (tunable parameter)

// coupling force on ball i (normalized by NUM_NEIGHBORS, so that its size does not grow with the coupling width):
//   F_i = 1/K * sum_{0 < |j-i| <= K} [ (c_i + c_j)/2 * phi(h_j - h_i) + COUPLING_DAMPING * (v_j - v_i) ]
//   phi(x) = x + SPRING_CUBIC * x^3 (hardening spring)
// impact of ball i (bounce or clamp, velocity change dv_i): v_j += IMPACT_TRANSFER / K * dv_i for its neighbors j
#define COUPLING_DAMPING 0.1
#define SPRING_CUBIC 1.0
#define IMPACT_TRANSFER 0.2
#endif

// Time events: every kick schedules the next one at
//   t_next = t + base_event_period * (1 + EVENT_HEIGHT_COUPLING * mean_i(pre(h_i)))
// so that the next event time depends on the state (LS-SA, event time sensitivities).
#define EVENT_HEIGHT_COUPLING 0.1
// lower bound for the period (relative to base_event_period), avoids event chattering
#define MIN_EVENT_PERIOD_FACTOR 1e-3

typedef uint32_t ValueReference;

// LS-SA: which branch of the event map a ball took in the last event iteration (the partial derivatives
// of the event map are evaluated lazily from this, see model.c)
#define EVENT_NONE    0   // not affected (identity)
#define EVENT_BOUNCE  1   // bounced: h+ = ground, v+ = -e * v-
#define EVENT_CLAMPED 2   // bounced and clamped: h+ = ground, v+ = 0

typedef struct {

    // continuous states
    double h;
    double v;

    // (tunable) parameters
    double e;
    double k;
    double drag_coefficient;

    // values right before the last event iteration, pre(h), pre(v)
    double pre_h;
    double pre_v;

    // discrete state
    int32_t bounce_count;

    // LS-SA: event branch taken in the last event iteration (EVENT_*), 1 byte instead of a cached Jacobian
    uint8_t eventKind;

#if NUM_NEIGHBORS > 0
    // spring stiffness (tunable parameter)
    double c;
#endif

} BallData; // 64 bytes (coupled: 72 bytes)

typedef struct {
    BallData balls[NUM_BALLS];

    double g;
    double v_min;
    double ground_flexibility;
    double base_event_period;

    double event_kick_velocity;

    // LS-SA: whether the last event iteration was a time event (kick), needed for the lazy event-map partials
    bool lastEventKick;

    // LS-SA: cached partials of nextEventTime, computed when the time event was handled
    // (all zero if the last event iteration did not (re-)schedule the next time event)
    double dTn_dtime; // d nextEventTime / d time (of the current event)
    double dTn_dpreh; // d nextEventTime / d pre(h_i), identical for all balls
    double dTn_dbase; // d nextEventTime / d base_event_period

#if NUM_NEIGHBORS > 0
    // balls that bounced / were clamped in the last event iteration (their impulse was passed to the neighbors);
    // LS-SA evaluates the coupling part of the event map from this list
    uint32_t impacted[NUM_BALLS];
    uint32_t nImpacted;
#endif

} ModelData;

#endif /* config_h */
