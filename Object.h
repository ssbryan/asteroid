#ifndef OBJECT_H
#define OBJECT_H

// Author Stephen Bryan
// June 22, 2017

#include <map>
#include <set>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <assert.h>

class ObjectMgr;
class Object;

const double G = 6.67408e-11;

extern double kMinVforDVCheck;
extern double kMaxDVChange;
extern double kMaxDVNormalChange;
extern double kStartDataSave;
extern unsigned int kNthDataSave;

struct LocFlags
{
    int mCoarse;
    int mMid;
    int mFine;
};

struct DeltaV
{
    DeltaV()
        : mDvx(0)
        , mDvy(0)
        , mDvz(0)
        , mObj(0)
    {}

    DeltaV(double dvx, double dvy, double dvz, Object* obj)
        : mDvx(dvx)
        , mDvy(dvy)
        , mDvz(dvz)
        , mObj(obj)
    {}

    DeltaV(const DeltaV& dv)
        : mDvx(dv.mDvx)
        , mDvy(dv.mDvy)
        , mDvz(dv.mDvz)
        , mObj(dv.mObj)
    {}

    double mDvx;
    double mDvy;
    double mDvz;
    Object* mObj;
};

struct OutputLoc
{
    OutputLoc(double x, double y, double z)
        : mX(x)
        , mY(y)
        , mZ(z)
    {}

    float mX;
    float mY;
    float mZ;
};

// data calculated per timestep
struct TimeStepData
{
    TimeStepData(double x, double y, double z, double vx, double vy, double vz)
        : mVx(vx)
        , mVy(vy)
        , mVz(vz)
    {
        mLoc[0] = x;
        mLoc[1] = y;
        mLoc[2] = z;
    }

    void    ClearData(void)
    {
        // clear all but deltaVs
        // we want to preserve them when we take mini timesteps
        mLoc[0] = 0;
        mLoc[1] = 0;
        mLoc[2] = 0;
        mVx = 0;
        mVy = 0;
        mVz = 0;
//        mDeltaV.clear();
//        mSubs.clear();
    }

    double  MagV2(void) const;
    double  MagDV2(double& dvx, double& dvy, double& dvz) const;

    // data
    // location, meters
    double          mLoc[3];

    // velocity, m/s
    double          mVx;
    double          mVy;
    double          mVz;

    std::map<int, DeltaV> mDeltaV;

    // constituents, if any
    std::map<int, Object*> mSubs;
};


class Object
{
public:
    Object(float x, float y, float z, double vx, double vy, double vz, double mass, int index);
    ~Object();

    void    GetLocation(double& x, double& y, double& z);
    void    GetVelocity(double& vx, double& vy, double& vz);
    void    SetLocation(double x, double y, double z);
    void    SetVelocity(double vx, double vy, double vz);
    double  GetMass() const;
    double  GetGMass() const;
    bool    IsDone(Object* obj, bool done);
    void    SetDone(Object* obj, bool done);
    int     GetIndex(void);
    double  GetRadius() const;
    void    CommitTimeStepData(void);
    void    ClearCurrentData(void);

    // write out position info for one timestep
    bool    WriteStep(const ObjectMgr* mgr, bool last);
    std::ofstream* GetOFile(void);
    std::string GetOFileName(void);
    void    AddDeltaV(double dvx, double dvy, double dvz, Object* obj);
    void    AppendDeltaV(double dvx, double dvy, double dvz, Object* obj);
    void    RemoveDeltaVandInitialize(Object* obj);
    void    ProcessDeltaVs(double timestep);
    void    ProcessFractionalDeltaVs(std::set<Object*>& multiGp, double timestep, double factor);
    bool    CalcLocFlags(double coarse, double xmin, double ymin, double zmin);
    bool    IntersectsInTimestep(Object* obj, 
                                 double& proximityTime,
                                 double& dist,
                                 double closeFactor,
                                 double currentTimestep);
    double  CheckTimestepCriteria(void) const;
    TimeStepData&   GetCommitted(void);
    const LocFlags& GetLocFlagsX(void) const;
    const LocFlags& GetLocFlagsY(void) const;
    const LocFlags& GetLocFlagsZ(void) const;
    const LocFlags& GetLocFlags2X(void) const;
    const LocFlags& GetLocFlags2Y(void) const;
    const LocFlags& GetLocFlags2Z(void) const;

private:
    // index for generating output filename
    int             mIndex;
    std::string     mFname;
    std::ofstream   mOfile;

    // mass * G, m^3/s^2
    double          mGMass;
    double          mMass;

    // deterined from mass
    double          mRadius;

    // changeable data
    TimeStepData    mCommitted;
    TimeStepData    mCurrent;
    double          mLastDeltaV;
    bool            mUseLastDV;

    // flags to determine what part of space the object is in
    LocFlags        mLocFlagsX;
    LocFlags        mLocFlagsY;
    LocFlags        mLocFlagsZ;
    // half shifted location
    LocFlags        mLocFlags2X;
    LocFlags        mLocFlags2Y;
    LocFlags        mLocFlags2Z;

    // flags for whether object has been evaluated for the mapped objects in this round
    std::map<Object*, bool> mDoneFlags;

    // committed data for output, in float form
    std::vector<OutputLoc> mOutput;
};

inline void Object::GetLocation(double& x, double& y, double& z)
{
    x = mCommitted.mLoc[0];
    y = mCommitted.mLoc[1];
    z = mCommitted.mLoc[2];
}

inline void Object::GetVelocity(double& vx, double& vy, double& vz)
{
    vx = mCommitted.mVx;
    vy = mCommitted.mVy;
    vz = mCommitted.mVz;
}

inline double Object::GetMass() const
{
    return mMass;
}

inline double Object::GetGMass() const
{
    return mGMass;
}

inline bool Object::IsDone(Object* obj, bool done)
{
    std::map<Object*, bool>::iterator oIt = mDoneFlags.find(obj);

    if (oIt != mDoneFlags.end())
    {
        return done == oIt->second;
    }

    return false;
}

inline void Object::SetDone(Object* obj, bool done)
{
    mDoneFlags[obj] = done;
}

inline void Object::SetLocation(double x, double y, double z)
{
    mCurrent.mLoc[0] = x;
    mCurrent.mLoc[1] = y;
    mCurrent.mLoc[2] = z;
}

inline void Object::SetVelocity(double vx, double vy, double vz)
{
    mCurrent.mVx = vx;
    mCurrent.mVy = vy;
    mCurrent.mVz = vz;
}

inline std::ofstream* Object::GetOFile(void)
{
    return &mOfile;
}

inline std::string Object::GetOFileName(void)
{
    return mFname;
}

inline int Object::GetIndex(void)
{
    return mIndex;
}

inline double Object::GetRadius() const
{
    return mRadius;
}

inline void Object::CommitTimeStepData()
{
    mCommitted = mCurrent;
}

inline void Object::ClearCurrentData()
{
    mCurrent.mDeltaV.clear();
    mCurrent.ClearData();
}

inline TimeStepData& Object::GetCommitted(void)
{
    return mCommitted;
}

inline const LocFlags& Object::GetLocFlagsX(void) const
{
    return mLocFlagsX;
}

inline const LocFlags& Object::GetLocFlagsY(void) const
{
    return mLocFlagsY;
}

inline const LocFlags& Object::GetLocFlagsZ(void) const
{
    return mLocFlagsZ;
}

inline const LocFlags& Object::GetLocFlags2X(void) const
{
    return mLocFlags2X;
}

inline const LocFlags& Object::GetLocFlags2Y(void) const
{
    return mLocFlags2Y;
}

inline const LocFlags& Object::GetLocFlags2Z(void) const
{
    return mLocFlags2Z;
}

#endif
