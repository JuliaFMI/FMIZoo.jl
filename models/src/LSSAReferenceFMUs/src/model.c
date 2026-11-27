#include <math.h>   // for fabs(), INFINITY
#include <float.h>  // for DBL_MIN
#include <stdlib.h> // for calloc(), realloc(), free()
#include <string.h> // for memcmp(), memcpy()
#include "model.h"

/*
 * N bouncing balls with LS-SA support.
 *
 * Per ball i (value references START_VR + (i-1)*VARS_PER_BALL + OFF_*):
 *   der(h) = v
 *   der(v) = g + k*time - drag*(1-e)*v*|v|
 *   event indicator  z = h - ground_level,  ground_level = -ground_flexibility*(1-e)
 *   error indicator    = -h
 *
 * Events (eventUpdate), in this order within one event iteration:
 *   pre(h) := h, pre(v) := v                      (for all balls)
 *   time event (time >= nextEventTime):           v := v + event_kick_velocity (all balls),
 *       nextEventTime := time + base_event_period * (1 + EVENT_HEIGHT_COUPLING * mean(pre(h)))
 *   state event (h <= ground_level && v < 0):     h := ground_level, v := -e*v (v := 0 if |v| < v_min)
 *
 * LS-SA:
 *   - fmi3GetDirectionalDerivative / fmi3GetAdjointDerivative are sparse, O(nUnknowns + nKnowns)
 *     (plus O(NUM_BALLS) for the inherently dense rows/columns of nextEventTime, time, g, ground_flexibility),
 *     with cached plans, streaming fast paths and struct-of-arrays scratch, see the LS-SA section below.
 *   - partials w.r.t. time (independent variable) are available in every mode
 *   - partials w.r.t. tunable parameters are available in every mode
 *   - event map: d x+ / d pre(x) is evaluated lazily from the branch each ball took in the last event iteration
 *     (event handling only records 1 byte per ball; identity if a ball was not affected);
 *     d x+ / d p (p = e, ground_flexibility) is returned in Event Mode only (in Continuous-Time Mode the
 *     outputs h, v do not depend on parameters)
 *   - d nextEventTime / d {pre(h_i), time, base_event_period} is cached when the time event is handled
 *     (zero if the last event iteration did not schedule a new time event)
 */

#define V_MIN (0.1)


/* ------------------------------------------------------------------------------------------------
 * model equations
 * --------------------------------------------------------------------------------------------- */

static inline double g_eff(ModelInstance* comp, const BallData* b) {
    return M(g) + b->k * comp->time;
}

static inline double derV(ModelInstance* comp, const BallData* b) {
    const double c = b->drag_coefficient * (1.0 - b->e);
    return g_eff(comp, b) - c * b->v * fabs(b->v);
}

static inline double groundLevel(ModelInstance* comp, const BallData* b) {
    return -M(ground_flexibility) * (1.0 - b->e);
}

double nextEventPeriod(ModelInstance* comp, double sum_h, bool* clamped) {
    const double P = M(base_event_period);
    const double period = P * (1.0 + EVENT_HEIGHT_COUPLING * sum_h / NUM_BALLS);
    const double minPeriod = P * MIN_EVENT_PERIOD_FACTOR;

    if (period < minPeriod) {
        *clamped = true;
        return minPeriod;
    }

    *clamped = false;
    return period;
}

// LS-SA: local Jacobian of the event map x+ = a(x-, p) of the last event iteration, evaluated lazily from
// the branch taken (b->eventKind) - event handling only records one byte per ball. (`comp` must be in scope.)
//   none:    h+ = h-,      v+ = v- (+ kick)
//   bounce:  h+ = ground,  v+ = -e * (v- (+ kick))      ground = -ground_flexibility * (1 - e)
//   clamped: h+ = ground,  v+ = 0
#define J_H_PREH(b) ((b)->eventKind == EVENT_NONE ? 1.0 : 0.0)                                   // d h+ / d pre(h)
#define J_V_PREV(b) ((b)->eventKind == EVENT_NONE ? 1.0 : ((b)->eventKind == EVENT_BOUNCE ? -(b)->e : 0.0)) // d v+ / d pre(v)
#define J_H_E(b)    ((b)->eventKind == EVENT_NONE ? 0.0 : M(ground_flexibility))                  // d h+ / d e
#define J_H_GF(b)   ((b)->eventKind == EVENT_NONE ? 0.0 : -(1.0 - (b)->e))                        // d h+ / d ground_flexibility
#define J_V_E(b)    ((b)->eventKind == EVENT_BOUNCE ? -((b)->pre_v + (M(lastEventKick) ? M(event_kick_velocity) : 0.0)) : 0.0) // d v+ / d e

void setStartValues(ModelInstance *comp) {
    M(g) = -9.81;
    M(v_min) = V_MIN;
    M(ground_flexibility) = 0.1;
    M(base_event_period) = 1.0;
    M(event_kick_velocity) = 2.0;

    M(lastEventKick) = false;
    M(dTn_dtime) = 0.0;
    M(dTn_dpreh) = 0.0;
    M(dTn_dbase) = 0.0;

    for (int i = 0; i < NUM_BALLS; i++) {
        BallData* b = &comp->modelData.balls[i];

        b->h = 10.0 + (double) i;
        b->v = 0.0;
        b->e = 0.7;
        b->k = 0.0;
        b->drag_coefficient = 0.1;
        b->bounce_count = 0;

        b->pre_h = b->h;
        b->pre_v = b->v;
        b->eventKind = EVENT_NONE;
    }
}

Status calculateValues(ModelInstance *comp) {
    UNUSED(comp);
    // nothing to do, all values are computed on demand
    return OK;
}

/* ------------------------------------------------------------------------------------------------
 * get / set
 * --------------------------------------------------------------------------------------------- */

Status getFloat64(ModelInstance* comp, ValueReference vr, double values[], size_t nValues, size_t* index) {

    ASSERT_NVALUES(1);

    if (vr < START_VR) {
        switch (vr) {
            case vr_time:
                values[(*index)++] = comp->time;
                return OK;
            case vr_g:
                values[(*index)++] = M(g);
                return OK;
            case vr_v_min:
                values[(*index)++] = M(v_min);
                return OK;
            case vr_ground_flexibility:
                values[(*index)++] = M(ground_flexibility);
                return OK;
            case vr_base_event_period:
                values[(*index)++] = M(base_event_period);
                return OK;
            case vr_nextEventTime:
                values[(*index)++] = comp->nextEventTimeDefined ? comp->nextEventTime : INFINITY;
                return OK;
            default:
                logError(comp, "Get Float64 is not allowed for value reference %u.", vr);
                return Error;
        }
    }

    const size_t offset = vr - START_VR;
    const size_t ballIdx = offset / VARS_PER_BALL;
    const size_t varType = offset % VARS_PER_BALL;

    if (ballIdx >= NUM_BALLS) {
        logError(comp, "Get Float64 is not allowed for value reference %u.", vr);
        return Error;
    }

    const BallData* b = &comp->modelData.balls[ballIdx];

    switch (varType) {
        case OFF_H:
            values[(*index)++] = b->h;
            return OK;
        case OFF_DER_H:
        case OFF_V:
            values[(*index)++] = b->v;
            return OK;
        case OFF_DER_V:
            values[(*index)++] = derV(comp, b);
            return OK;
        case OFF_E:
            values[(*index)++] = b->e;
            return OK;
        case OFF_K:
            values[(*index)++] = b->k;
            return OK;
        case OFF_DRAG:
            values[(*index)++] = b->drag_coefficient;
            return OK;
        case OFF_PRE_H:
            values[(*index)++] = b->pre_h;
            return OK;
        case OFF_PRE_V:
            values[(*index)++] = b->pre_v;
            return OK;
        case OFF_ERR_H:
            values[(*index)++] = -(b->h);
            return OK;
        case OFF_EVENT_G:
            values[(*index)++] = b->h - groundLevel(comp, b);
            return OK;
        default:
            logError(comp, "Get Float64 is not allowed for value reference %u.", vr);
            return Error;
    }
}

// FMI 3.0: tunable parameters may only be changed in Event Mode (ME) / Step Mode (CS), not in Continuous-Time Mode
static bool tunableSettable(ModelInstance* comp, ValueReference vr) {
#if FMI_VERSION == 3
    if (comp->state & (Instantiated | InitializationMode | ConfigurationMode | ReconfigurationMode | EventMode | StepMode)) {
        return true;
    }
    logError(comp, "Tunable parameter with value reference %u can only be set in Event Mode (or Step Mode), not in Continuous-Time Mode.", vr);
    return false;
#else
    UNUSED(comp);
    UNUSED(vr);
    return true;
#endif
}

Status setFloat64(ModelInstance* comp, ValueReference vr, double value) {

    if (vr < START_VR) {
        switch (vr) {
            case vr_g:
#if FMI_VERSION > 1
                if (comp->type == ModelExchange &&
                    comp->state != Instantiated &&
                    comp->state != InitializationMode) {
                    logError(comp, "Variable g can only be set after instantiation or in initialization mode.");
                    return Error;
                }
#endif
                M(g) = value;
                return OK;

            case vr_ground_flexibility:
                if (!tunableSettable(comp, vr)) return Error;
                M(ground_flexibility) = value;
                return OK;

            case vr_base_event_period:
                if (!tunableSettable(comp, vr)) return Error;
                M(base_event_period) = value;
                return OK;

            case vr_v_min:
                 M(v_min) = value;
                 return OK;

            default:
                logError(comp, "Cannot set global VR %u", vr);
                return Error;
        }
    }

    const size_t offset = vr - START_VR;
    const size_t ballIdx = offset / VARS_PER_BALL;
    const size_t varType = offset % VARS_PER_BALL;

    if (ballIdx >= NUM_BALLS) {
        logError(comp, "Unexpected value reference: %u.", vr);
        return Error;
    }

    BallData* b = &comp->modelData.balls[ballIdx];

    switch (varType) {
        case OFF_H:
            b->h = value;
            return OK;

        case OFF_V:
            b->v = value;
            return OK;

        case OFF_E:
            if (!tunableSettable(comp, vr)) return Error;
            b->e = value;
            return OK;

        case OFF_K:
            if (!tunableSettable(comp, vr)) return Error;
            b->k = value;
            return OK;

        case OFF_DRAG:
            if (!tunableSettable(comp, vr)) return Error;
            b->drag_coefficient = value;
            return OK;

        default:
            logError(comp, "Unexpected value reference: %u.", vr);
            return Error;
    }
}

Status getOutputDerivative(ModelInstance *comp, ValueReference valueReference, int order, double *value) {

    if (order != 1) {
        logError(comp, "The output derivative order %d for value reference %u is not available.", order, valueReference);
        return Error;
    }

    if (valueReference < START_VR) {
        logError(comp, "The output derivative for value reference %u is not available.", valueReference);
        return Error;
    }

    const size_t offset = valueReference - START_VR;
    const size_t ballIdx = offset / VARS_PER_BALL;
    const size_t varType = offset % VARS_PER_BALL;

    if (ballIdx >= NUM_BALLS) return Error;
    const BallData* b = &comp->modelData.balls[ballIdx];

    switch (varType) {
    case OFF_H:
        *value = b->v;
        return OK;
    case OFF_V:
        *value = derV(comp, b);
        return OK;
    default:
        logError(comp, "The output derivative for value reference %u is not available.", valueReference);
        return Error;
    }
}

Status getInt32(ModelInstance* comp, ValueReference vr,
                int32_t values[], size_t nValues, size_t* index) {
    ASSERT_NVALUES(1);

    if (vr < START_VR) return Error;

    const size_t offset = vr - START_VR;
    const size_t ballIdx = offset / VARS_PER_BALL;
    const size_t varType = offset % VARS_PER_BALL;

    if (ballIdx >= NUM_BALLS) return Error;

    if (varType == OFF_BOUNCE) {
        values[(*index)++] = comp->modelData.balls[ballIdx].bounce_count;
        return OK;
    }

    logError(comp, "Get Int32 is not allowed for value reference %u.", vr);
    return Error;
}

Status setInt32(ModelInstance* comp, ValueReference vr, const int32_t values[], size_t nValues, size_t* index)
{
    ASSERT_NVALUES(1);
    if (vr < START_VR) return Error;

    const size_t offset = vr - START_VR;
    const size_t ballIdx = offset / VARS_PER_BALL;
    const size_t varType = offset % VARS_PER_BALL;

    if (ballIdx >= NUM_BALLS) return Error;

    if (varType == OFF_BOUNCE) {
        // LS-SA: discrete states are settable in Event Mode
        if (comp->state != EventMode) {
            logError(comp, "Discrete state with value reference %u can only be set in Event Mode.", vr);
            return Error;
        }
        comp->modelData.balls[ballIdx].bounce_count = values[(*index)++];
        return OK;
    }

    logError(comp, "Set Int32 is not allowed for value reference %u.", vr);
    return Error;
}

/* ------------------------------------------------------------------------------------------------
 * continuous part
 * --------------------------------------------------------------------------------------------- */

size_t getNumberOfEventIndicators(ModelInstance* comp) {
    UNUSED(comp);
    return NUM_BALLS;
}

size_t getNumberOfContinuousStates(ModelInstance* comp) {
    UNUSED(comp);
    return NUM_BALLS * 2;
}

Status getContinuousStates(ModelInstance *comp, double x[], size_t nx) {
    UNUSED(nx);
    for (int i = 0; i < NUM_BALLS; i++) {
        const BallData* b = &comp->modelData.balls[i];
        x[2 * i]     = b->h;
        x[2 * i + 1] = b->v;
    }
    return OK;
}

Status setContinuousStates(ModelInstance *comp, const double x[], size_t nx) {
    UNUSED(nx);
    for (int i = 0; i < NUM_BALLS; i++) {
        BallData* b = &comp->modelData.balls[i];
        b->h = x[2 * i];
        b->v = x[2 * i + 1];
    }
    return OK;
}

Status getDerivatives(ModelInstance *comp, double dx[], size_t nx) {
    UNUSED(nx);
    if (!dx) return OK;
    for (int i = 0; i < NUM_BALLS; i++) {
        const BallData* b = &comp->modelData.balls[i];
        dx[2 * i]     = b->v;
        dx[2 * i + 1] = derV(comp, b);
    }
    return OK;
}

Status getEventIndicators(ModelInstance *comp, double z[], size_t nz) {
    UNUSED(nz);
    for (int i = 0; i < NUM_BALLS; i++) {
        const BallData* b = &comp->modelData.balls[i];
        z[i] = b->h - groundLevel(comp, b);
    }
    return OK;
}

/* ------------------------------------------------------------------------------------------------
 * events
 * --------------------------------------------------------------------------------------------- */

Status eventUpdate(ModelInstance *comp) {

    comp->valuesOfContinuousStatesChanged = false;
    comp->nominalsOfContinuousStatesChanged = false;
    comp->terminateSimulation = false;
    comp->newDiscreteStatesNeeded = false;

    const bool timeEvent = comp->nextEventTimeDefined && comp->time >= comp->nextEventTime;
    const double gf = M(ground_flexibility);
    double sum_pre_h = 0.0;

    for (int i = 0; i < NUM_BALLS; i++) {
        BallData* b = &comp->modelData.balls[i];

        // values right before this event iteration
        b->pre_h = b->h;
        b->pre_v = b->v;
        sum_pre_h += b->h;

        // LS-SA: only the branch taken is recorded, the event-map partials are evaluated lazily (J_* below)
        b->eventKind = EVENT_NONE;

        // time event: kick
        if (timeEvent) {
            b->v += M(event_kick_velocity);
            comp->valuesOfContinuousStatesChanged = true;
        }

        // state event: bounce
        const double ground = -gf * (1.0 - b->e);

        if (b->h <= ground && b->v < 0) {
            b->bounce_count++;
            b->h = ground + DBL_MIN;
            b->v = -b->e * b->v;
            b->eventKind = EVENT_BOUNCE;

            if (fabs(b->v) < M(v_min)) {
                b->v = 0.0;
                b->eventKind = EVENT_CLAMPED;
            }

            comp->valuesOfContinuousStatesChanged = true;
        }
    }

    M(lastEventKick) = timeEvent;

    if (timeEvent) {
        if (M(base_event_period) > 0) {
            bool clamped = false;
            const double period = nextEventPeriod(comp, sum_pre_h, &clamped);

            comp->nextEventTime = comp->time + period;

            M(dTn_dtime) = 1.0;
            if (clamped) {
                M(dTn_dpreh) = 0.0;
                M(dTn_dbase) = MIN_EVENT_PERIOD_FACTOR;
            } else {
                M(dTn_dpreh) = M(base_event_period) * EVENT_HEIGHT_COUPLING / NUM_BALLS;
                M(dTn_dbase) = 1.0 + EVENT_HEIGHT_COUPLING * sum_pre_h / NUM_BALLS;
            }
        } else {
            comp->nextEventTimeDefined = false;
            M(dTn_dtime) = 0.0;
            M(dTn_dpreh) = 0.0;
            M(dTn_dbase) = 0.0;
        }
    } else {
        // nextEventTime was not (re-)computed in this event iteration
        M(dTn_dtime) = 0.0;
        M(dTn_dpreh) = 0.0;
        M(dTn_dbase) = 0.0;
    }

    return OK;
}

/* ------------------------------------------------------------------------------------------------
 * LS-SA: directional (JVP) and adjoint (VJP) derivatives
 *
 * FMI passes explicit value reference lists with every call. To keep the per-call overhead close to
 * the unavoidable cost of reading these lists once, three standard techniques are used:
 *
 *  1. Plans: every (unknowns, knowns) pair is analysed once (validation + pattern detection) and
 *     cached per instance (PLAN_CACHE_SIZE entries, keyed by the list contents via memcmp). Importers
 *     call the same pattern in every step, only the seed changes. Build with LSSA_PLAN_CACHE=OFF to
 *     analyse every call instead (used to measure "cold" calls).
 *  2. Fast paths: lists that are block structured over the balls (e.g. h1,v1,h2,v2,..., all event
 *     indicators, all parameters e_i) or consist of a single global variable are processed by
 *     streaming kernels reading the seeds directly with fixed strides: no scatter / clear, no per
 *     element value reference decoding.
 *  3. Struct of arrays: for arbitrary lists (generic path) the seeds are scattered into one contiguous
 *     array per variable type (cache friendly), the outputs are evaluated, the slots are cleared again.
 *
 * All paths use the same per-ball formulas (jvpLocal / vjpLocal, transposes of each other), so they
 * cannot diverge. Every variable depends on itself (identity); everything else is the analytic partial
 * derivative of the model equations or a cached event partial.
 * --------------------------------------------------------------------------------------------- */

#ifndef PLAN_CACHE_SIZE
#define PLAN_CACHE_SIZE 8
#endif

// generic path scratch: globals, then one contiguous array of NUM_BALLS entries per variable type
#define SCRATCH_SIZE (NUM_GLOBAL_VRS + (size_t)VARS_PER_BALL * NUM_BALLS)

typedef enum { LIST_GENERIC, LIST_BLOCK, LIST_GLOBAL } ListKind;

typedef struct {
    ListKind kind;
    size_t m;                     // BLOCK: variables per ball (group size)
    size_t ball0;                 // BLOCK: first ball (0-based)
    size_t nBalls;                // BLOCK: number of consecutive balls
    int types[VARS_PER_BALL];     // BLOCK: variable type (OFF_*) of slot s within a group
    ValueReference vr;            // GLOBAL: the value reference
} ListPattern;

typedef struct {
    bool used;
    size_t nUnknowns;
    size_t nKnowns;
    size_t capacity;
    ValueReference* vrs;          // copy of the unknowns followed by the knowns
    ListPattern unknowns;
    ListPattern knowns;
} Plan;

typedef struct {
    double* scratch;              // all zero between calls
    Plan plans[PLAN_CACHE_SIZE];
    size_t nextPlan;              // round robin replacement
} LSSAData;

Status allocateModelBuffers(ModelInstance* comp) {
    LSSAData* d = (LSSAData*)calloc(1, sizeof(LSSAData));
    if (!d) return Error;
    d->scratch = (double*)calloc(SCRATCH_SIZE, sizeof(double));
    if (!d->scratch) {
        free(d);
        return Error;
    }
    comp->lssa = d;
    return OK;
}

void freeModelBuffers(ModelInstance* comp) {
    LSSAData* d = (LSSAData*)comp->lssa;
    if (!d) return;
    for (size_t i = 0; i < PLAN_CACHE_SIZE; i++) {
        free(d->plans[i].vrs);
    }
    free(d->scratch);
    free(d);
    comp->lssa = NULL;
}

static inline size_t scratchIndex(ValueReference vr) {
    if (vr < START_VR) return (size_t)vr;
    const size_t offset = vr - START_VR;
    return NUM_GLOBAL_VRS + (offset % VARS_PER_BALL) * NUM_BALLS + offset / VARS_PER_BALL;
}

// all Float64 variables (bounce_count is Int32 and therefore not differentiable)
static bool isFloat64VR(ValueReference vr) {
    if (vr < NUM_GLOBAL_VRS) return true;
    if (vr < START_VR) return false;
    const size_t offset = vr - START_VR;
    if (offset >= (size_t)NUM_BALLS * VARS_PER_BALL) return false;
    return (offset % VARS_PER_BALL) != OFF_BOUNCE;
}

static Status checkVRs(ModelInstance* comp, const char* f, const char* arg, const ValueReference vrs[], size_t n) {
    if (n > 0 && !vrs) {
        logError(comp, "%s: Argument %s must not be NULL.", f, arg);
        return Error;
    }
    for (size_t i = 0; i < n; i++) {
        if (!isFloat64VR(vrs[i])) {
            logError(comp, "%s: %s[%zu] = %u is not a Float64 variable of this FMU.", f, arg, i, vrs[i]);
            return Error;
        }
    }
    return OK;
}

/* ---- per-ball formulas ------------------------------------------------------------------------ */

#define SB(t) sb[(size_t)(t) * st]

// d (unknown of variable type `type` of ball b) * seed; seeds of the knowns of this ball at sb[t * st],
// seeds of the global knowns at sg[vr]
static inline double jvpLocal(ModelInstance* comp, const BallData* b, size_t type,
                              const double* sb, size_t st, const double* sg, bool eventMode) {
    switch (type) {
        case OFF_H: {
            double d = SB(OFF_H) + J_H_PREH(b) * SB(OFF_PRE_H);
            if (eventMode) {
                d += J_H_E(b) * SB(OFF_E) + J_H_GF(b) * sg[vr_ground_flexibility];
            }
            return d;
        }
        case OFF_V: {
            double d = SB(OFF_V) + J_V_PREV(b) * SB(OFF_PRE_V);
            if (eventMode) {
                d += J_V_E(b) * SB(OFF_E);
            }
            return d;
        }
        case OFF_DER_H:
            return SB(OFF_DER_H) + SB(OFF_V);
        case OFF_DER_V: {
            const double c   = b->drag_coefficient * (1.0 - b->e);
            const double vav = b->v * fabs(b->v);
            return SB(OFF_DER_V)
                 - 2.0 * c * fabs(b->v)      * SB(OFF_V)
                 + b->drag_coefficient * vav * SB(OFF_E)
                 - (1.0 - b->e) * vav        * SB(OFF_DRAG)
                 + comp->time                * SB(OFF_K)
                 + sg[vr_g]
                 + b->k                      * sg[vr_time];
        }
        case OFF_ERR_H:
            return SB(OFF_ERR_H) - SB(OFF_H);
        case OFF_EVENT_G:
            return SB(OFF_EVENT_G) + SB(OFF_H)
                 - M(ground_flexibility) * SB(OFF_E)
                 + (1.0 - b->e) * sg[vr_ground_flexibility];
        default:
            return SB(type); // e, k, drag, pre(h), pre(v): identity
    }
}

#undef SB
#define WB(t) wb[(size_t)(t) * st]

// seed * d unknowns / d (known of variable type `type` of ball b); seeds of the unknowns of this ball at
// wb[t * st], seeds of the global unknowns at wg[vr]
static inline double vjpLocal(ModelInstance* comp, const BallData* b, size_t type,
                              const double* wb, size_t st, const double* wg, bool eventMode) {
    switch (type) {
        case OFF_H:
            return WB(OFF_H) - WB(OFF_ERR_H) + WB(OFF_EVENT_G);
        case OFF_V: {
            const double c = b->drag_coefficient * (1.0 - b->e);
            return WB(OFF_V) + WB(OFF_DER_H) - 2.0 * c * fabs(b->v) * WB(OFF_DER_V);
        }
        case OFF_E: {
            double d = WB(OFF_E)
                     + b->drag_coefficient * b->v * fabs(b->v) * WB(OFF_DER_V)
                     - M(ground_flexibility) * WB(OFF_EVENT_G);
            if (eventMode) {
                d += J_H_E(b) * WB(OFF_H) + J_V_E(b) * WB(OFF_V);
            }
            return d;
        }
        case OFF_K:
            return WB(OFF_K) + comp->time * WB(OFF_DER_V);
        case OFF_DRAG:
            return WB(OFF_DRAG) - (1.0 - b->e) * b->v * fabs(b->v) * WB(OFF_DER_V);
        case OFF_PRE_H:
            return WB(OFF_PRE_H) + J_H_PREH(b) * WB(OFF_H) + M(dTn_dpreh) * wg[vr_nextEventTime];
        case OFF_PRE_V:
            return WB(OFF_PRE_V) + J_V_PREV(b) * WB(OFF_V);
        default:
            return WB(type); // der(h), der(v), error / event indicator: identity
    }
}

// contribution of the unknowns of ball b to the adjoint of the global known vr
static inline double vjpGlobalFromBall(ModelInstance* comp, const BallData* b, ValueReference vr,
                                       const double* wb, size_t st, bool eventMode) {
    UNUSED(comp);
    switch (vr) {
        case vr_time:
            return b->k * WB(OFF_DER_V);
        case vr_g:
            return WB(OFF_DER_V);
        case vr_ground_flexibility: {
            double d = (1.0 - b->e) * WB(OFF_EVENT_G);
            if (eventMode) {
                d += J_H_GF(b) * WB(OFF_H);
            }
            return d;
        }
        default:
            return 0.0;
    }
}

#undef WB

static inline bool globalHasBallTerms(ValueReference vr) {
    return vr == vr_time || vr == vr_g || vr == vr_ground_flexibility;
}

// global unknown: identity + dependence on global knowns
static inline double jvpGlobalSelf(ModelInstance* comp, ValueReference vr, const double* sg) {
    double d = sg[vr];
    if (vr == vr_nextEventTime) {
        d += M(dTn_dtime) * sg[vr_time] + M(dTn_dbase) * sg[vr_base_event_period];
    }
    return d;
}

// global known: identity + contribution of global unknowns
static inline double vjpGlobalSelf(ModelInstance* comp, ValueReference vr, const double* wg) {
    double d = wg[vr];
    if (vr == vr_time) {
        d += M(dTn_dtime) * wg[vr_nextEventTime];
    } else if (vr == vr_base_event_period) {
        d += M(dTn_dbase) * wg[vr_nextEventTime];
    }
    return d;
}

/* ---- plans ------------------------------------------------------------------------------------ */

// block structure: groups of m consecutive entries per ball, the same variables for consecutive balls
// (validates implicitly: matching lists only contain Float64 variables of existing balls)
static bool detectBlock(const ValueReference vrs[], size_t n, ListPattern* p) {

    if (n == 0 || vrs[0] < START_VR) return false;

    // first group: consecutive entries of the same ball
    const size_t ball0 = (vrs[0] - START_VR) / VARS_PER_BALL;
    size_t m = 0;
    while (m < n && m < VARS_PER_BALL && vrs[m] >= START_VR && (vrs[m] - START_VR) / VARS_PER_BALL == ball0) {
        m++;
    }
    if (n % m != 0) return false;

    const size_t nBalls = n / m;
    if (ball0 + nBalls > NUM_BALLS) return false;

    bool seen[VARS_PER_BALL] = { false };
    for (size_t s = 0; s < m; s++) {
        const size_t t = (vrs[s] - START_VR) % VARS_PER_BALL;
        if (seen[t] || t == OFF_BOUNCE) return false; // duplicates / non-Float64: generic path (validates)
        seen[t] = true;
        p->types[s] = (int)t;
    }

    // all further groups: the same variables of the following balls
    for (size_t g = 1; g < nBalls; g++) {
        const ValueReference* grp = vrs + g * m;
        const ValueReference shift = (ValueReference)(g * VARS_PER_BALL);
        for (size_t s = 0; s < m; s++) {
            if (grp[s] != vrs[s] + shift) return false;
        }
    }

    p->m = m;
    p->ball0 = ball0;
    p->nBalls = nBalls;
    return true;
}

static Status analyzeList(ModelInstance* comp, const char* f, const char* arg,
                          const ValueReference vrs[], size_t n, ListPattern* p) {

    if (n > 0 && !vrs) {
        logError(comp, "%s: Argument %s must not be NULL.", f, arg);
        return Error;
    }

    if (n == 1 && vrs[0] < NUM_GLOBAL_VRS) {
        p->kind = LIST_GLOBAL;
        p->vr = vrs[0];
        return OK;
    }

    if (detectBlock(vrs, n, p)) {
        p->kind = LIST_BLOCK;
        return OK;
    }

    p->kind = LIST_GENERIC;
    return checkVRs(comp, f, arg, vrs, n);
}

static Status getPlan(ModelInstance* comp, const char* f,
                      const ValueReference unknowns[], size_t nUnknowns,
                      const ValueReference knowns[], size_t nKnowns,
                      Plan* tmp, const Plan** plan) {
#ifdef LSSA_PLAN_CACHE
    LSSAData* d = (LSSAData*)comp->lssa;
    for (size_t i = 0; i < PLAN_CACHE_SIZE; i++) {
        const Plan* p = &d->plans[i];
        if (p->used && p->nUnknowns == nUnknowns && p->nKnowns == nKnowns &&
            (nUnknowns == 0 || memcmp(p->vrs, unknowns, nUnknowns * sizeof(ValueReference)) == 0) &&
            (nKnowns == 0 || memcmp(p->vrs + nUnknowns, knowns, nKnowns * sizeof(ValueReference)) == 0)) {
            *plan = p;
            return OK;
        }
    }
#endif

    if (analyzeList(comp, f, "unknowns", unknowns, nUnknowns, &tmp->unknowns) != OK) return Error;
    if (analyzeList(comp, f, "knowns", knowns, nKnowns, &tmp->knowns) != OK) return Error;
    tmp->nUnknowns = nUnknowns;
    tmp->nKnowns = nKnowns;
    *plan = tmp;

#ifdef LSSA_PLAN_CACHE
    const size_t n = nUnknowns + nKnowns;
    Plan* p = &d->plans[d->nextPlan];
    if (n > p->capacity) {
        ValueReference* mem = (ValueReference*)realloc(p->vrs, n * sizeof(ValueReference));
        if (!mem) return OK; // not cached, still correct
        p->vrs = mem;
        p->capacity = n;
    }
    if (nUnknowns > 0) memcpy(p->vrs, unknowns, nUnknowns * sizeof(ValueReference));
    if (nKnowns > 0) memcpy(p->vrs + nUnknowns, knowns, nKnowns * sizeof(ValueReference));
    p->nUnknowns = nUnknowns;
    p->nKnowns = nKnowns;
    p->unknowns = tmp->unknowns;
    p->knowns = tmp->knowns;
    p->used = true;
    d->nextPlan = (d->nextPlan + 1) % PLAN_CACHE_SIZE;
    *plan = p;
#endif

    return OK;
}

/* ---- generic path (arbitrary lists, struct-of-arrays scratch) ---------------------------------- */

static void jvpGeneric(ModelInstance* comp, const ValueReference unknowns[], size_t nUnknowns,
                       const ValueReference knowns[], size_t nKnowns,
                       const double seed[], double sensitivity[], bool eventMode) {
    double* s = ((LSSAData*)comp->lssa)->scratch;

    for (size_t j = 0; j < nKnowns; j++) {
        s[scratchIndex(knowns[j])] += seed[j];
    }

    for (size_t i = 0; i < nUnknowns; i++) {
        const ValueReference vr = unknowns[i];
        if (vr < START_VR) {
            double d = jvpGlobalSelf(comp, vr, s);
            if (vr == vr_nextEventTime && M(dTn_dpreh) != 0.0) {
                const double* sp = s + NUM_GLOBAL_VRS + (size_t)OFF_PRE_H * NUM_BALLS;
                double acc = 0.0;
                for (size_t k = 0; k < NUM_BALLS; k++) {
                    acc += sp[k];
                }
                d += M(dTn_dpreh) * acc;
            }
            sensitivity[i] = d;
        } else {
            const size_t offset = vr - START_VR;
            const size_t ball = offset / VARS_PER_BALL;
            sensitivity[i] = jvpLocal(comp, &comp->modelData.balls[ball], offset % VARS_PER_BALL,
                                      s + NUM_GLOBAL_VRS + ball, NUM_BALLS, s, eventMode);
        }
    }

    for (size_t j = 0; j < nKnowns; j++) {
        s[scratchIndex(knowns[j])] = 0.0;
    }
}

static void vjpGeneric(ModelInstance* comp, const ValueReference unknowns[], size_t nUnknowns,
                       const ValueReference knowns[], size_t nKnowns,
                       const double seed[], double sensitivity[], bool eventMode) {
    double* w = ((LSSAData*)comp->lssa)->scratch;

    for (size_t i = 0; i < nUnknowns; i++) {
        w[scratchIndex(unknowns[i])] += seed[i];
    }

    for (size_t j = 0; j < nKnowns; j++) {
        const ValueReference vr = knowns[j];
        if (vr < START_VR) {
            double d = vjpGlobalSelf(comp, vr, w);
            if (globalHasBallTerms(vr)) {
                for (size_t k = 0; k < NUM_BALLS; k++) {
                    d += vjpGlobalFromBall(comp, &comp->modelData.balls[k], vr, w + NUM_GLOBAL_VRS + k, NUM_BALLS, eventMode);
                }
            }
            sensitivity[j] = d;
        } else {
            const size_t offset = vr - START_VR;
            const size_t ball = offset / VARS_PER_BALL;
            sensitivity[j] = vjpLocal(comp, &comp->modelData.balls[ball], offset % VARS_PER_BALL,
                                      w + NUM_GLOBAL_VRS + ball, NUM_BALLS, w, eventMode);
        }
    }

    for (size_t i = 0; i < nUnknowns; i++) {
        w[scratchIndex(unknowns[i])] = 0.0;
    }
}

/* ---- fast paths (block structured lists): type-major streaming kernels -------------------------
 * One tight loop over the balls per variable type of the output list (the switch is outside the loop),
 * the seeds are read in place through per-type column pointers / strides (absent types: a zero with
 * stride 0). Ball ranges of unknowns and knowns may differ: the output range is split into at most
 * three segments (without / with / without seeds), so the loops contain no membership tests.
 * The formulas are the same as in jvpLocal / vjpLocal (checked by the tests: fast == generic path). */

static const double ZERO = 0.0;

typedef struct {
    const double* col[VARS_PER_BALL];  // per variable type: seed of ball `ball0`
    size_t stride[VARS_PER_BALL];      // per variable type: distance between consecutive balls
    size_t ball0;
} Columns;

static void columnsZero(Columns* c) {
    for (size_t t = 0; t < VARS_PER_BALL; t++) {
        c->col[t] = &ZERO;
        c->stride[t] = 0;
    }
    c->ball0 = 0;
}

static void columnsFromBlock(const ListPattern* p, const double seed[], Columns* c) {
    columnsZero(c);
    for (size_t s = 0; s < p->m; s++) {
        c->col[p->types[s]] = seed + s;
        c->stride[p->types[s]] = p->m;
    }
    c->ball0 = p->ball0;
}

// overlap [lo, hi) of the ball ranges [a0, a1) and [b0, b1); empty overlap: lo = hi = a1
static void overlap(size_t a0, size_t a1, size_t b0, size_t b1, size_t* lo, size_t* hi) {
    *lo = a0 > b0 ? a0 : b0;
    *hi = a1 < b1 ? a1 : b1;
    if (*lo >= *hi) {
        *lo = a1;
        *hi = a1;
    }
}

// column access (pointer and stride hoisted into locals, so they cannot alias with the output)
#define COLUMN(name, t) const double* const name##_p = c->col[t]; const size_t name##_s = c->stride[t]
#define AT(name) name##_p[(i - b0) * name##_s]

// unknowns of variable type `type` of balls [i0, i1) -> out[(i - i0) * os]; seeds of the knowns in c
static void jvpKernel(ModelInstance* comp, size_t type, size_t i0, size_t i1, const Columns* c,
                      const double* sg, bool eventMode, double* out, size_t os) {
    const BallData* balls = comp->modelData.balls;
    const size_t b0 = c->ball0;

    switch (type) {
        case OFF_H: {
            COLUMN(h, OFF_H); COLUMN(ph, OFF_PRE_H); COLUMN(e, OFF_E);
            const double sgf = sg[vr_ground_flexibility];
            for (size_t i = i0; i < i1; i++, out += os) {
                const BallData* b = &balls[i];
                double d = AT(h) + J_H_PREH(b) * AT(ph);
                if (eventMode) d += J_H_E(b) * AT(e) + J_H_GF(b) * sgf;
                *out = d;
            }
            break;
        }
        case OFF_V: {
            COLUMN(v, OFF_V); COLUMN(pv, OFF_PRE_V); COLUMN(e, OFF_E);
            for (size_t i = i0; i < i1; i++, out += os) {
                const BallData* b = &balls[i];
                double d = AT(v) + J_V_PREV(b) * AT(pv);
                if (eventMode) d += J_V_E(b) * AT(e);
                *out = d;
            }
            break;
        }
        case OFF_DER_H: {
            COLUMN(dh, OFF_DER_H); COLUMN(v, OFF_V);
            for (size_t i = i0; i < i1; i++, out += os) {
                *out = AT(dh) + AT(v);
            }
            break;
        }
        case OFF_DER_V: {
            COLUMN(dv, OFF_DER_V); COLUMN(v, OFF_V); COLUMN(e, OFF_E); COLUMN(dr, OFF_DRAG); COLUMN(k, OFF_K);
            const double sgg = sg[vr_g], sgt = sg[vr_time], t = comp->time;
            for (size_t i = i0; i < i1; i++, out += os) {
                const BallData* b = &balls[i];
                const double av = fabs(b->v), vav = b->v * av;
                *out = AT(dv)
                     - 2.0 * b->drag_coefficient * (1.0 - b->e) * av * AT(v)
                     + b->drag_coefficient * vav * AT(e)
                     - (1.0 - b->e) * vav * AT(dr)
                     + t * AT(k)
                     + sgg
                     + b->k * sgt;
            }
            break;
        }
        case OFF_ERR_H: {
            COLUMN(er, OFF_ERR_H); COLUMN(h, OFF_H);
            for (size_t i = i0; i < i1; i++, out += os) {
                *out = AT(er) - AT(h);
            }
            break;
        }
        case OFF_EVENT_G: {
            COLUMN(z, OFF_EVENT_G); COLUMN(h, OFF_H); COLUMN(e, OFF_E);
            const double gf = M(ground_flexibility), sgf = sg[vr_ground_flexibility];
            for (size_t i = i0; i < i1; i++, out += os) {
                *out = AT(z) + AT(h) - gf * AT(e) + (1.0 - balls[i].e) * sgf;
            }
            break;
        }
        default: { // e, k, drag, pre(h), pre(v): identity
            COLUMN(x, type);
            for (size_t i = i0; i < i1; i++, out += os) {
                *out = AT(x);
            }
            break;
        }
    }
}

// knowns of variable type `type` of balls [i0, i1) -> out[(i - i0) * os]; seeds of the unknowns in c,
// seeds of the global unknowns in wg
static void vjpKernel(ModelInstance* comp, size_t type, size_t i0, size_t i1, const Columns* c,
                      const double* wg, bool eventMode, double* out, size_t os) {
    const BallData* balls = comp->modelData.balls;
    const size_t b0 = c->ball0;

    switch (type) {
        case OFF_H: {
            COLUMN(h, OFF_H); COLUMN(er, OFF_ERR_H); COLUMN(z, OFF_EVENT_G);
            for (size_t i = i0; i < i1; i++, out += os) {
                *out = AT(h) - AT(er) + AT(z);
            }
            break;
        }
        case OFF_V: {
            COLUMN(v, OFF_V); COLUMN(dh, OFF_DER_H); COLUMN(dv, OFF_DER_V);
            for (size_t i = i0; i < i1; i++, out += os) {
                const BallData* b = &balls[i];
                *out = AT(v) + AT(dh) - 2.0 * b->drag_coefficient * (1.0 - b->e) * fabs(b->v) * AT(dv);
            }
            break;
        }
        case OFF_E: {
            COLUMN(e, OFF_E); COLUMN(dv, OFF_DER_V); COLUMN(z, OFF_EVENT_G); COLUMN(h, OFF_H); COLUMN(v, OFF_V);
            const double gf = M(ground_flexibility);
            for (size_t i = i0; i < i1; i++, out += os) {
                const BallData* b = &balls[i];
                double d = AT(e) + b->drag_coefficient * b->v * fabs(b->v) * AT(dv) - gf * AT(z);
                if (eventMode) d += J_H_E(b) * AT(h) + J_V_E(b) * AT(v);
                *out = d;
            }
            break;
        }
        case OFF_K: {
            COLUMN(k, OFF_K); COLUMN(dv, OFF_DER_V);
            const double t = comp->time;
            for (size_t i = i0; i < i1; i++, out += os) {
                *out = AT(k) + t * AT(dv);
            }
            break;
        }
        case OFF_DRAG: {
            COLUMN(dr, OFF_DRAG); COLUMN(dv, OFF_DER_V);
            for (size_t i = i0; i < i1; i++, out += os) {
                const BallData* b = &balls[i];
                *out = AT(dr) - (1.0 - b->e) * b->v * fabs(b->v) * AT(dv);
            }
            break;
        }
        case OFF_PRE_H: {
            COLUMN(ph, OFF_PRE_H); COLUMN(h, OFF_H);
            const double tn = M(dTn_dpreh) * wg[vr_nextEventTime];
            for (size_t i = i0; i < i1; i++, out += os) {
                *out = AT(ph) + J_H_PREH(&balls[i]) * AT(h) + tn;
            }
            break;
        }
        case OFF_PRE_V: {
            COLUMN(pv, OFF_PRE_V); COLUMN(v, OFF_V);
            for (size_t i = i0; i < i1; i++, out += os) {
                *out = AT(pv) + J_V_PREV(&balls[i]) * AT(v);
            }
            break;
        }
        default: { // der(h), der(v), error / event indicator: identity
            COLUMN(x, type);
            for (size_t i = i0; i < i1; i++, out += os) {
                *out = AT(x);
            }
            break;
        }
    }
}

// contributions of the unknowns of balls [i0, i1) to the adjoint of the global known vr
static double vjpGlobalKernel(ModelInstance* comp, ValueReference vr, size_t i0, size_t i1,
                              const Columns* c, bool eventMode) {
    const BallData* balls = comp->modelData.balls;
    const size_t b0 = c->ball0;
    double d = 0.0;

    switch (vr) {
        case vr_time: {
            COLUMN(dv, OFF_DER_V);
            for (size_t i = i0; i < i1; i++) d += balls[i].k * AT(dv);
            break;
        }
        case vr_g: {
            COLUMN(dv, OFF_DER_V);
            for (size_t i = i0; i < i1; i++) d += AT(dv);
            break;
        }
        case vr_ground_flexibility: {
            COLUMN(z, OFF_EVENT_G); COLUMN(h, OFF_H);
            for (size_t i = i0; i < i1; i++) {
                d += (1.0 - balls[i].e) * AT(z);
                if (eventMode) d += J_H_GF(&balls[i]) * AT(h);
            }
            break;
        }
        default:
            break;
    }
    return d;
}

#undef AT
#undef COLUMN

static bool jvpFast(ModelInstance* comp, const Plan* plan, const double seed[], double sensitivity[], bool eventMode) {
    const ListPattern* U = &plan->unknowns;
    const ListPattern* K = &plan->knowns;
    double sg[NUM_GLOBAL_VRS] = { 0.0 };
    Columns none, cols;
    columnsZero(&none);
    cols = none;

    if (U->kind == LIST_BLOCK && (K->kind == LIST_BLOCK || K->kind == LIST_GLOBAL)) {
        const size_t u0 = U->ball0, u1 = U->ball0 + U->nBalls;
        size_t lo = u1, hi = u1;  // balls of the unknowns with seeded per-ball knowns
        if (K->kind == LIST_BLOCK) {
            columnsFromBlock(K, seed, &cols);
            overlap(u0, u1, K->ball0, K->ball0 + K->nBalls, &lo, &hi);
        } else {
            sg[K->vr] = seed[0];
        }
        for (size_t s = 0; s < U->m; s++) {
            const size_t t = (size_t)U->types[s];
            double* out = sensitivity + s;
            jvpKernel(comp, t, u0, lo, &none, sg, eventMode, out, U->m);
            jvpKernel(comp, t, lo, hi, &cols, sg, eventMode, out + (lo - u0) * U->m, U->m);
            jvpKernel(comp, t, hi, u1, &none, sg, eventMode, out + (hi - u0) * U->m, U->m);
        }
        return true;
    }

    if (U->kind == LIST_GLOBAL && K->kind == LIST_BLOCK) {
        // global unknowns do not depend on per-ball knowns, except nextEventTime on pre(h)
        double d = 0.0;
        if (U->vr == vr_nextEventTime && M(dTn_dpreh) != 0.0) {
            for (size_t s = 0; s < K->m; s++) {
                if (K->types[s] == OFF_PRE_H) {
                    const double* sp = seed + s;
                    double acc = 0.0;
                    for (size_t g = 0; g < K->nBalls; g++) {
                        acc += sp[g * K->m];
                    }
                    d = M(dTn_dpreh) * acc;
                }
            }
        }
        sensitivity[0] = d;
        return true;
    }

    return false;
}

static bool vjpFast(ModelInstance* comp, const Plan* plan, const double seed[], double sensitivity[], bool eventMode) {
    const ListPattern* U = &plan->unknowns;
    const ListPattern* K = &plan->knowns;
    double wg[NUM_GLOBAL_VRS] = { 0.0 };
    Columns none, cols;
    columnsZero(&none);
    cols = none;

    if (K->kind == LIST_BLOCK && (U->kind == LIST_BLOCK || U->kind == LIST_GLOBAL)) {
        const size_t k0 = K->ball0, k1 = K->ball0 + K->nBalls;
        size_t lo = k1, hi = k1;  // balls of the knowns with seeded per-ball unknowns
        if (U->kind == LIST_BLOCK) {
            columnsFromBlock(U, seed, &cols);
            overlap(k0, k1, U->ball0, U->ball0 + U->nBalls, &lo, &hi);
        } else {
            wg[U->vr] = seed[0];
        }
        for (size_t s = 0; s < K->m; s++) {
            const size_t t = (size_t)K->types[s];
            double* out = sensitivity + s;
            vjpKernel(comp, t, k0, lo, &none, wg, eventMode, out, K->m);
            vjpKernel(comp, t, lo, hi, &cols, wg, eventMode, out + (lo - k0) * K->m, K->m);
            vjpKernel(comp, t, hi, k1, &none, wg, eventMode, out + (hi - k0) * K->m, K->m);
        }
        return true;
    }

    if (K->kind == LIST_GLOBAL && U->kind == LIST_BLOCK) {
        // the global known is not among the per-ball unknowns: only their contributions
        columnsFromBlock(U, seed, &cols);
        sensitivity[0] = vjpGlobalKernel(comp, K->vr, U->ball0, U->ball0 + U->nBalls, &cols, eventMode);
        return true;
    }

    return false;
}

/* ---- entry points ------------------------------------------------------------------------------ */

Status getDirectionalDerivative(ModelInstance* comp,
    const ValueReference unknowns[], size_t nUnknowns,
    const ValueReference knowns[], size_t nKnowns,
    const double seed[], size_t nSeed,
    double sensitivity[], size_t nSensitivity) {

    static const char* f = "fmi3GetDirectionalDerivative";

    if (nSeed != nKnowns) {
        logError(comp, "%s: nSeed = %zu must equal nKnowns = %zu.", f, nSeed, nKnowns);
        return Error;
    }
    if (nSensitivity != nUnknowns) {
        logError(comp, "%s: nSensitivity = %zu must equal nUnknowns = %zu.", f, nSensitivity, nUnknowns);
        return Error;
    }
    if ((nKnowns > 0 && !seed) || (nUnknowns > 0 && !sensitivity)) {
        logError(comp, "%s: Arguments seed / sensitivity must not be NULL.", f);
        return Error;
    }

    Plan tmp;
    const Plan* plan = NULL;
    if (getPlan(comp, f, unknowns, nUnknowns, knowns, nKnowns, &tmp, &plan) != OK) return Error;

    const bool eventMode = (comp->state == EventMode);
    if (!jvpFast(comp, plan, seed, sensitivity, eventMode)) {
        jvpGeneric(comp, unknowns, nUnknowns, knowns, nKnowns, seed, sensitivity, eventMode);
    }
    return OK;
}

Status getAdjointDerivative(ModelInstance* comp,
    const ValueReference unknowns[], size_t nUnknowns,
    const ValueReference knowns[], size_t nKnowns,
    const double seed[], size_t nSeed,
    double sensitivity[], size_t nSensitivity) {

    static const char* f = "fmi3GetAdjointDerivative";

    if (nSeed != nUnknowns) {
        logError(comp, "%s: nSeed = %zu must equal nUnknowns = %zu.", f, nSeed, nUnknowns);
        return Error;
    }
    if (nSensitivity != nKnowns) {
        logError(comp, "%s: nSensitivity = %zu must equal nKnowns = %zu.", f, nSensitivity, nKnowns);
        return Error;
    }
    if ((nUnknowns > 0 && !seed) || (nKnowns > 0 && !sensitivity)) {
        logError(comp, "%s: Arguments seed / sensitivity must not be NULL.", f);
        return Error;
    }

    Plan tmp;
    const Plan* plan = NULL;
    if (getPlan(comp, f, unknowns, nUnknowns, knowns, nKnowns, &tmp, &plan) != OK) return Error;

    const bool eventMode = (comp->state == EventMode);
    if (!vjpFast(comp, plan, seed, sensitivity, eventMode)) {
        vjpGeneric(comp, unknowns, nUnknowns, knowns, nKnowns, seed, sensitivity, eventMode);
    }
    return OK;
}

// single Jacobian entry (used by the FMI 2.0 interface and the dense fallback)
Status getPartialDerivative(ModelInstance *comp, ValueReference unknown, ValueReference known, double *pd) {

    if (!isFloat64VR(unknown) || !isFloat64VR(known)) {
        logError(comp, "Partial derivative of %u w.r.t. %u is not defined (not a Float64 variable).", unknown, known);
        return Error;
    }

    const double seed = 1.0;
    jvpGeneric(comp, &unknown, 1, &known, 1, &seed, pd, comp->state == EventMode);
    return OK;
}
