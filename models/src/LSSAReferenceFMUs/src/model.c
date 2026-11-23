#include <math.h>   // for fabs()
#include <float.h>  // for DBL_MIN
#include <string.h> // for memset()
#include "model.h"


#ifndef NUM_BALLS
#define NUM_BALLS 1
#endif

#define START_VR 100
#define VARS_PER_BALL 12

#define V_MIN (0.1)
#define EVENT_EPSILON (1e-10)

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

#define VR_TIME 0
#define VR_G 1
#define VR_V_MIN 2
#define VR_GROUND_FLEX 3
#define VR_BASE_PERIOD 4
#define VR_NEXT_EVENT 5

static inline double g_eff(ModelInstance* comp, BallData* b)
{
    return M(g) + (b->k * comp->time);
}

static inline void checkGlobalDirty(ModelInstance *comp) {
    if (comp->isDirtyJacobian) {
        comp->jacobianEpoch++;
        comp->isDirtyJacobian = false;
    }
}

void markJacobianDirty(ModelInstance *comp) {
    comp->jacobianEpoch++;
    comp->isDirtyJacobian = false;
}

static inline void updateBallTimeCache(ModelInstance *comp, size_t ballIdx) {
    checkGlobalDirty(comp);
    BallJacobian* J = &comp->ballJacobians[ballIdx];
    if (J->epochTime == comp->jacobianEpoch) return;

    BallData* b = &comp->modelData.balls[ballIdx];
    double v = b->v;
    double c = b->drag_coefficient * (1.0 - b->e);

    J->dh_dt     = v;
    J->dv_dt     = (M(g) + b->k * comp->time) - (c * v * fabs(v));
    J->dderv_dt  = b->k;
    J->epochTime = comp->jacobianEpoch;
}

void setStartValues(ModelInstance *comp) {
    M(g) = -9.81;
    M(v_min) = V_MIN;
    M(ground_flexibility) = 0.1;
    M(base_event_period) = 1.0;
    M(current_event_period) = M(base_event_period) * 0.7;
    M(event_kick_velocity) = 2.0;
    M(nextEventTime) = M(base_event_period);
    comp->nextEventTime = M(base_event_period);

    comp->jacobianEpoch = 1;
    comp->isDirtyJacobian = false;

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
        b->der_h = 0;
        b->der_v = 0;

        BallJacobian* J = &comp->ballJacobians[i];
        J->epochTime   = 0;
        J->epochDerV   = 0;
        J->epochEventG = 0;

        // Preload fixed matrix entries once at startup:
        // 1) Event indicator z with respect to state x = (h, v)
        J->dz_dh = 1.0;
        J->dz_dv = 0.0;

        // 2) State derivative x. = (der_h, der_v) with respect to state x = (h, v)
        J->dderh_dh = 0.0;
        J->dderh_dv = 1.0;
        J->dderv_dh = 0.0;
    }
}

Status calculateValues(ModelInstance *comp) {
    UNUSED(comp);
    // nothing to do
    return OK;
}

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
                logError(comp, "Global GetFloat64 unknown VR %u.", vr);
                return Error;
        }
    }

    size_t offset = vr - START_VR;
    size_t ballIdx = offset / VARS_PER_BALL;
    size_t varType = offset % VARS_PER_BALL;

    if (ballIdx >= NUM_BALLS) return Error;

    BallData* b = &comp->modelData.balls[ballIdx];

    switch (varType) {
        case OFF_H:
            values[(*index)++] = b->h;
            return OK;
        case OFF_DER_H:
        case OFF_V:
            values[(*index)++] = b->v;
            return OK;
        case OFF_DER_V:
            updateBallTimeCache(comp, ballIdx);
            values[(*index)++] = comp->ballJacobians[ballIdx].dv_dt;
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
            values[(*index)++] = b->h - (-(M(ground_flexibility) * (1.0 - b->e)));
            return OK;
        default:
            logError(comp, "Get Float64 is not allowed for value reference %u.", vr);
            return Error;
    }
}


Status setFloat64(ModelInstance* comp, ValueReference vr, double value) {

    if (vr < START_VR) {
        comp->jacobianEpoch++;
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
                M(ground_flexibility) = value;
                return OK;

            case vr_base_event_period:
                M(base_event_period) = value;
                M(current_event_period) = value * 0.7;
                return OK;
                
            case vr_v_min:
                 M(v_min) = value;
                 return OK;

            default: 
                logError(comp, "Cannot set global VR %u", vr);
                return Error;
        }
    }

    size_t offset = vr - START_VR;
    size_t ballIdx = offset / VARS_PER_BALL;
    size_t varType = offset % VARS_PER_BALL;

    if (ballIdx >= NUM_BALLS) return Error;
    
    BallData* b = &comp->modelData.balls[ballIdx];
    BallJacobian* J = &comp->ballJacobians[ballIdx];

    // Mark only the affected ball as dirty in O(1) (Epoch = 0)!
    switch (varType) {
        case OFF_H:
            b->h = value;
            return OK;

        case OFF_V:
            b->v = value;
            J->epochTime = 0;
            J->epochDerV = 0;
            return OK;

        case OFF_E:
            b->e = value;
            J->epochTime   = 0;
            J->epochDerV   = 0;
            J->epochEventG = 0;
            return OK;

        case OFF_K:
            b->k = value;
            J->epochTime = 0;
            return OK;

        case OFF_DRAG:
            b->drag_coefficient = value;
            J->epochTime = 0;
            J->epochDerV = 0;
            return OK;

        case OFF_PRE_H:
            b->pre_h = value;
            return OK;

        case OFF_PRE_V:
            b->pre_v = value;
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

    if (valueReference < START_VR) return Error;

    size_t offset = valueReference - START_VR;
    size_t ballIdx = offset / VARS_PER_BALL;
    size_t varType = offset % VARS_PER_BALL;

    if (ballIdx >= NUM_BALLS) return Error;
    BallData* b = &comp->modelData.balls[ballIdx];

    switch (varType) {
    case OFF_H:
        *value = b->v;
        return OK;
    case OFF_V:
        updateBallTimeCache(comp, ballIdx);
        *value = comp->ballJacobians[ballIdx].dv_dt;
        return OK;
    default:
        logError(comp, "The output derivative for value reference %u is not available.", valueReference);
        return Error;
    }
}

void updateBallDerVCache(ModelInstance *comp, size_t ballIdx) {
    checkGlobalDirty(comp);
    BallJacobian* J = &comp->ballJacobians[ballIdx];
    if (J->epochDerV == comp->jacobianEpoch) return;

    BallData* b = &comp->modelData.balls[ballIdx];
    double v = b->v;
    double abs_v = fabs(v);
    double v_abs_v = v * abs_v;
    double one_minus_e = 1.0 - b->e;
    double c = b->drag_coefficient * one_minus_e;

    J->dderv_dk    = comp->time;
    J->dderv_dv    = -2.0 * c * abs_v;
    J->dderv_de    = b->drag_coefficient * v_abs_v;
    J->dderv_ddrag = -one_minus_e * v_abs_v;
    J->epochDerV   = comp->jacobianEpoch;
}

static inline void updateBallEventGCache(ModelInstance *comp, size_t ballIdx) {
    checkGlobalDirty(comp);
    BallJacobian* J = &comp->ballJacobians[ballIdx];
    if (J->epochEventG == comp->jacobianEpoch) return;

    BallData* b = &comp->modelData.balls[ballIdx];
    J->devg_de     = -M(ground_flexibility);
    J->devg_dflex  = 1.0 - b->e;
    J->epochEventG = comp->jacobianEpoch;
}

// Special cases (event indicator z with respect to state x; state derivative x. with respect to state x) using preloaded matrix

// 1. Event indicator z (N) with respect to state x (2N)
static inline bool isEventIndicatorWrtState(const ValueReference unknowns[], size_t nUnknowns, const ValueReference knowns[], size_t nKnowns) {
    if (nUnknowns != NUM_BALLS || nKnowns != 2 * NUM_BALLS) return false;
    return (unknowns[0] == START_VR + OFF_EVENT_G) &&
           (knowns[0]   == START_VR + OFF_H) &&
           (knowns[1]   == START_VR + OFF_V);
}

// 2. State derivative x. (2N) with respect to state x (2N)
static inline bool isDerStateWrtState(const ValueReference unknowns[], size_t nUnknowns, const ValueReference knowns[], size_t nKnowns) {
    if (nUnknowns != 2 * NUM_BALLS || nKnowns != 2 * NUM_BALLS) return false;
    return (unknowns[0] == START_VR + OFF_DER_H) &&
           (unknowns[1] == START_VR + OFF_DER_V) &&
           (knowns[0]   == START_VR + OFF_H) &&
           (knowns[1]   == START_VR + OFF_V);
}

bool evalSpecialJVP(ModelInstance *comp, const ValueReference unknowns[], size_t nUnknowns, const ValueReference knowns[], size_t nKnowns, const double seed[], double sensitivity[]) {
    // Special case 1: Event indicator z with respect to state x
    // Fixed matrix (dz_dh = 1.0, dz_dv = 0.0) is already preloaded in setStartValues -> no cache update needed!
    if (isEventIndicatorWrtState(unknowns, nUnknowns, knowns, nKnowns)) {
        for (size_t i = 0; i < NUM_BALLS; i++) {
            const BallJacobian* J = &comp->ballJacobians[i];
            sensitivity[i] = J->dz_dh * seed[2 * i] + J->dz_dv * seed[2 * i + 1];
        }
        return true;
    }

    // Special case 2: State derivative x. with respect to state x
    // Fixed matrix entries (dderh_dh = 0.0, dderh_dv = 1.0, dderv_dh = 0.0) are preloaded;
    // only dderv_dv is lazily recomputed if the ball is dirty!
    if (isDerStateWrtState(unknowns, nUnknowns, knowns, nKnowns)) {
        checkGlobalDirty(comp);
        for (size_t i = 0; i < NUM_BALLS; i++) {
            BallJacobian* J = &comp->ballJacobians[i];
            const double s_h = seed[2 * i];
            const double s_v = seed[2 * i + 1];

            if (s_v != 0.0 && J->epochDerV != comp->jacobianEpoch) {
                updateBallDerVCache(comp, i);
            }

            sensitivity[2 * i]     = J->dderh_dh * s_h + J->dderh_dv * s_v;
            sensitivity[2 * i + 1] = J->dderv_dh * s_h + (s_v != 0.0 ? J->dderv_dv * s_v : 0.0);
        }
        return true;
    }

    return false;
}

bool evalSpecialVJP(ModelInstance *comp, const ValueReference unknowns[], size_t nUnknowns, const ValueReference knowns[], size_t nKnowns, const double seed[], double sensitivity[]) {
    // Special case 1: Event indicator z with respect to state x (transposed)
    // Fixed matrix (dz_dh = 1.0, dz_dv = 0.0) is already preloaded in setStartValues -> no cache update needed!
    if (isEventIndicatorWrtState(unknowns, nUnknowns, knowns, nKnowns)) {
        for (size_t i = 0; i < NUM_BALLS; i++) {
            const BallJacobian* J = &comp->ballJacobians[i];
            const double s_z = seed[i];
            sensitivity[2 * i]     = J->dz_dh * s_z;
            sensitivity[2 * i + 1] = J->dz_dv * s_z;
        }
        return true;
    }

    // Special case 2: State derivative x. with respect to state x (transposed)
    // Fixed matrix entries (dderh_dh = 0.0, dderh_dv = 1.0, dderv_dh = 0.0) are preloaded;
    // only dderv_dv is lazily recomputed if the ball is dirty!
    if (isDerStateWrtState(unknowns, nUnknowns, knowns, nKnowns)) {
        checkGlobalDirty(comp);
        for (size_t i = 0; i < NUM_BALLS; i++) {
            BallJacobian* J = &comp->ballJacobians[i];
            const double s_derh = seed[2 * i];
            const double s_derv = seed[2 * i + 1];

            if (s_derv != 0.0 && J->epochDerV != comp->jacobianEpoch) {
                updateBallDerVCache(comp, i);
            }

            sensitivity[2 * i]     = J->dderh_dh * s_derh + J->dderv_dh * s_derv;
            sensitivity[2 * i + 1] = J->dderh_dv * s_derh + (s_derv != 0.0 ? J->dderv_dv * s_derv : 0.0);
        }
        return true;
    }

    return false;
}

// Do not update directly; only update the requested cache (only update if dirty and needed!)
Status getPartialDerivative(ModelInstance *comp, ValueReference unknown, ValueReference known, double *pd)
{
    bool isUnkGlobal = (unknown < START_VR);
    bool isKnownGlobal = (known < START_VR);

    size_t idxUnk = 0, typeUnk = 0;
    if (!isUnkGlobal) {
        size_t offset = unknown - START_VR;
        idxUnk = offset / VARS_PER_BALL;
        typeUnk = offset % VARS_PER_BALL;
        if (idxUnk >= NUM_BALLS) return Error;
    }

    size_t idxKnown = 0, typeKnown = 0;
    if (!isKnownGlobal) {
        size_t offset = known - START_VR;
        idxKnown = offset / VARS_PER_BALL;
        typeKnown = offset % VARS_PER_BALL;
        if (idxKnown >= NUM_BALLS) return Error;
    }

    // Different balls do not affect each other -> immediately return 0!
    if (!isUnkGlobal && !isKnownGlobal && (idxUnk != idxKnown)) {
        *pd = 0.0;
        return OK;
    }

    // Case 1: Derivative with respect to time (known == vr_time)
    if (known == vr_time) {
        if (isUnkGlobal) {
            switch (unknown) {
                case vr_g:             *pd = 0.0; return OK;
                case vr_nextEventTime: *pd = 1.0; return OK;
                default:               *pd = 0.0; return OK;
            }
        } else {
            switch (typeUnk) {
                case OFF_H:
                    updateBallTimeCache(comp, idxUnk);
                    *pd = comp->ballJacobians[idxUnk].dh_dt;
                    return OK;
                case OFF_V:
                case OFF_DER_H:
                    updateBallTimeCache(comp, idxUnk);
                    *pd = comp->ballJacobians[idxUnk].dv_dt;
                    return OK;
                case OFF_DER_V:
                    updateBallTimeCache(comp, idxUnk);
                    *pd = comp->ballJacobians[idxUnk].dderv_dt;
                    return OK;
                default:
                    *pd = 0.0;
                    return OK;
            }
        }
    }

    // Case 2: Local Unknowns (Ball variables)
    if (!isUnkGlobal) {
        BallJacobian* J = &comp->ballJacobians[idxUnk];
        switch (typeUnk) {
            case OFF_H:
                *pd = (!isKnownGlobal && typeKnown == OFF_H) ? 1.0 : 0.0;
                return OK;
            case OFF_V:
                *pd = (!isKnownGlobal && typeKnown == OFF_V) ? 1.0 : 0.0;
                return OK;
            case OFF_DER_H:
                if (!isKnownGlobal) {
                    if (typeKnown == OFF_H) { *pd = J->dderh_dh; return OK; }
                    if (typeKnown == OFF_V) { *pd = J->dderh_dv; return OK; }
                }
                *pd = 0.0;
                return OK;
            case OFF_DER_V:
                if (isKnownGlobal) {
                    *pd = (known == vr_g) ? 1.0 : 0.0;
                    return OK;
                }
                switch (typeKnown) {
                    case OFF_H:
                        *pd = J->dderv_dh;
                        return OK;
                    case OFF_K:
                        updateBallDerVCache(comp, idxUnk);
                        *pd = J->dderv_dk;
                        return OK;
                    case OFF_V:
                        updateBallDerVCache(comp, idxUnk);
                        *pd = J->dderv_dv;
                        return OK;
                    case OFF_E:
                        updateBallDerVCache(comp, idxUnk);
                        *pd = J->dderv_de;
                        return OK;
                    case OFF_DRAG:
                        updateBallDerVCache(comp, idxUnk);
                        *pd = J->dderv_ddrag;
                        return OK;
                    default:
                        *pd = 0.0;
                        return OK;
                }
            case OFF_EVENT_G:
                if (isKnownGlobal) {
                    if (known == vr_ground_flexibility) {
                        updateBallEventGCache(comp, idxUnk);
                        *pd = J->devg_dflex;
                        return OK;
                    }
                    *pd = 0.0;
                    return OK;
                }
                if (typeKnown == OFF_H) { *pd = J->dz_dh; return OK; }
                if (typeKnown == OFF_V) { *pd = J->dz_dv; return OK; }
                if (typeKnown == OFF_E) {
                    updateBallEventGCache(comp, idxUnk);
                    *pd = J->devg_de;
                    return OK;
                }
                *pd = 0.0;
                return OK;
            case OFF_PRE_H:
            case OFF_PRE_V:
                *pd = 0.0;
                return OK;
            case OFF_ERR_H:
                *pd = (!isKnownGlobal && typeKnown == OFF_H) ? -1.0 : 0.0;
                return OK;
            default:
                *pd = 0.0;
                return OK;
        }
    }

    // Case 3: Global Unknowns
    if (unknown == vr_nextEventTime && known == vr_base_event_period) {
        *pd = 1.0;
        return OK;
    }

    *pd = 0.0;
    return OK;
}

Status eventUpdate(ModelInstance *comp) {
    comp->valuesOfContinuousStatesChanged = false;
    comp->nominalsOfContinuousStatesChanged = false;
    comp->terminateSimulation = false;

    if (comp->nextEventTimeDefined && comp->time >= comp->nextEventTime)
    {
        for (int i = 0; i < NUM_BALLS; i++) {
            BallData* b = &comp->modelData.balls[i];
            b->pre_v = b->v;
            b->v += M(event_kick_velocity);
        }

        comp->nextEventTime += M(current_event_period);
        comp->valuesOfContinuousStatesChanged = true;
    }

    for (int i = 0; i < NUM_BALLS; i++) {
        BallData* b = &comp->modelData.balls[i];

        double ground_level = -M(ground_flexibility) * (1.0 - b->e);

        if (b->h <= ground_level && b->v < 0) {
            b->bounce_count++;

            b->pre_h = b->h;
            b->pre_v = b->v;

            b->h = ground_level + DBL_MIN;

            b->v = -b->v * b->e;

            if (fabs(b->v) < M(v_min)) {
                b->v = 0.0;
            }

            comp->valuesOfContinuousStatesChanged = true;
        }
    }

    comp->jacobianEpoch++;

    return OK;
}

size_t getNumberOfEventIndicators(ModelInstance* comp) {

    UNUSED(comp);

    return NUM_BALLS;
}

size_t getNumberOfContinuousStates(ModelInstance* comp) {

    UNUSED(comp);

    return  NUM_BALLS * 2;
}

Status getContinuousStates(ModelInstance *comp, double x[], size_t nx) {

    UNUSED(nx);

    for (int i = 0; i < NUM_BALLS; i++) {
            BallData* b = &comp->modelData.balls[i];
            x[2 * i]     = b->h;  // Even indices = position
            x[2 * i + 1] = b->v;  // Uneven indices = speed
        }
        return OK;
}

Status setContinuousStates(ModelInstance *comp, const double x[], size_t nx) {

    UNUSED(nx);

    for (int i = 0; i < NUM_BALLS; i++) {
        BallData* b = &comp->modelData.balls[i];
        b->h = x[2 * i];      // update position
        b->v = x[2 * i + 1];  // update speed
    }
    comp->jacobianEpoch++;
    return OK;
}

Status getDerivatives(ModelInstance *comp, double dx[], size_t nx) {

    UNUSED(nx);

    for (int i = 0; i < NUM_BALLS; i++) {
        BallData* b = &comp->modelData.balls[i];
        b->der_h = b->v;

        double current_g = g_eff(comp, b);

        double effective_drag = b->drag_coefficient * (1.0 - b->e);
        double drag_accel = effective_drag * b->v * fabs(b->v);

        b->der_v = current_g - drag_accel;

        if (dx) {
            dx[2 * i]     = b->der_h;
            dx[2 * i + 1] = b->der_v;
        }
    }
    return OK;
}

Status getEventIndicators(ModelInstance *comp, double z[], size_t nz) {

    UNUSED(nz);

    for (int i = 0; i < NUM_BALLS; i++) {
        BallData* b = &comp->modelData.balls[i];
        double ground_level = -M(ground_flexibility) * (1.0 - b->e);
        z[i] = b->h - ground_level;
    }

    return OK;
}

Status getInt32(ModelInstance* comp, ValueReference vr,
                int32_t values[], size_t nValues, size_t* index) {
    ASSERT_NVALUES(1);
    
    if (vr < START_VR) return Error;

    size_t offset = vr - START_VR;
    size_t ballIdx = offset / VARS_PER_BALL;
    size_t varType = offset % VARS_PER_BALL;

    if (ballIdx >= NUM_BALLS) return Error;

    if (varType == OFF_BOUNCE) {
        values[(*index)++] = comp->modelData.balls[ballIdx].bounce_count;
        return OK;
    }
    
    return Error;
}

Status setInt32(ModelInstance* comp, ValueReference vr, const int32_t values[], size_t nValues, size_t* index)
{
    ASSERT_NVALUES(1);
    if (vr < START_VR) return Error;

    size_t offset = vr - START_VR;
    size_t ballIdx = offset / VARS_PER_BALL;
    size_t varType = offset % VARS_PER_BALL;

    if (ballIdx >= NUM_BALLS) return Error;

    if (varType == OFF_BOUNCE) {
         if (comp->state != EventMode) return Error;
         comp->modelData.balls[ballIdx].bounce_count = values[(*index)++];
         return OK;
    }
    return Error;
}