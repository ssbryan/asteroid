#ifndef OBJECTMGR_H
#define OBJECTMGR_H

// Author Stephen Bryan
// June 22, 2017

#include "Object.h"
#include <vector>
#include <string>

struct CollisionData
{
    CollisionData(Object* obj0, Object* obj1)
        : mObj0(obj0)
        , mObj1(obj1)
        , mTimeToCollision(-1)
        , mDist(0)
    {}

    CollisionData(const CollisionData& cd)
        : mObj0(cd.mObj0)
        , mObj1(cd.mObj1)
        , mTimeToCollision(cd.mTimeToCollision)
        , mDist(cd.mDist)
    {
    }

    bool    ContainsObject(Object* obj) const
    {
        if ((obj == mObj0) || (obj == mObj1))
        {
            return true;
        }

        return false;
    }

    double  mDist;
    double  mTimeToCollision;
    Object* mObj0;
    Object* mObj1;
};


class ObjectMgr
{
public:
    ObjectMgr(const std::string& fname);
    ~ObjectMgr(void);

    bool    Initialize(const std::string& fname);
    void    SetTimeData(double ts, double tend);
    bool    Run(void);
    bool    CalcStep(bool done);
    void    CalcStepForMultiObjects(std::vector<CollisionData>& multicolls, double tstep);
    void    CalcCollision(std::vector<CollisionData>& multicolls, double tstep, double tsfactor);
    double  GetOrbitTime(void) const;
    double  GetStepTime(void) const;
    void    GetRunTimeData(double& tstop, double& tend) const;
    void    ProcessCollisions(std::vector<CollisionData>& collisions);
    void    LogSimpleMessage(const char* msg);

private:
    std::ofstream   mLogFile;
    bool    mOkay;
    double  mTStep;
    double  mTEnd;
    double  mTStop;
    double  mPrint;
    double  mMaxOrbitTime;
    std::vector<Object*> mObjs;

    // extent
    double  mXmax;
    double  mXmin;
    double  mYmax;
    double  mYmin;
    double  mZmax;
    double  mZmin;

    // time
    time_t  mStartTime;
};

inline ObjectMgr::~ObjectMgr(void)
{
    for (std::vector<Object*>::iterator oIt = mObjs.begin();
        oIt != mObjs.end();
        ++oIt)
    {
        delete *oIt;
    }
}

inline void ObjectMgr::SetTimeData(double ts, double tend)
{
    mTStep = ts;
    mTEnd = tend;
    mTStop = mTEnd;
}

inline double ObjectMgr::GetOrbitTime(void) const
{
    return mMaxOrbitTime;
}

inline void ObjectMgr::GetRunTimeData(double& tstop, double& tend) const
{
    tstop = mTStop;
    tend = mTEnd;
}

inline double ObjectMgr::GetStepTime(void) const
{
    return mTStep;
}

#endif
