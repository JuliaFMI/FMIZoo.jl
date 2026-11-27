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
#define VARS_PER_BALL 12

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

} BallData; // 64 bytes

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

} ModelData;

#endif /* config_h */
