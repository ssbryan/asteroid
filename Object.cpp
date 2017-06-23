
// to get M_PI
#define _USE_MATH_DEFINES 

#include "Object.h"
#include "ObjectMgr.h"
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <assert.h>

#include <cmath>

// Author Stephen Bryan
// June 22, 2017

// tolerance values set in Options line
double kMinVforDVCheck = 0.01;
double kMaxDVChange = 0.01;
double kMaxDVNormalChange = 0.01;
double kStartDataSave = 0;
unsigned int kNthDataSave = 1;

const int maxskip = 10;


Object::Object(float x, float y, float z, double vx, double vy, double vz, double mass, int index)
    : mMass(mass)
    , mGMass(mass * G)
    , mFname("obj")
    , mIndex(index)
    , mCommitted(x, y, z, vx, vy, vz)
    , mCurrent(0, 0, 0, 0, 0, 0)
    , mLastDeltaV(0)
    , mUseLastDV(false)
{
    // determine radius: rho = M / V, V = M / rho

    // rho = 3g/cc plus porosity discount of 50%, plus random element that increases 
    double rho0 = 3000; // 3g/cc (* 1e-3kg/g * 1e6cc/m^3), similar to carbonaceous chondrites
    double porosity = 0.5; // 50% porosity, changes effective rho

    // should do an "if" to see if a random element is needed (solid rock, not rubble), 
    // or if we want a nickel-iron chunk at 7.9g/cc (3.84% of meteorites)
    double rho = rho0 * porosity;

    // volume
    double v = mMass / rho;

    // R = (3 * V / (4 * pi))^(1/3)
    mRadius = pow(3 * v / (4 * M_PI), 0.33333333) ;

#ifndef PC
    char num[9];
    sprintf(num, "%d", index);

    mFname += num;
    mFname += ".pos";
#else
    mFname += std::to_string(index) + ".pos";
#endif
    // clean out any existing same-named file
#ifndef PC
    mOfile.open(mFname.c_str(), std::ios_base::trunc | std::ios_base::out);
#else
    mOfile.open(mFname, std::ios::trunc | std::ios::out);
#endif
    mOfile.close();

    // open the file, read in initialization data and create objects
    // format of file (all doubles):  x y z\n
#ifndef PC
    mOfile.open(mFname.c_str(), std::ios_base::app | std::ios_base::out);
#else
    mOfile.open(mFname, std::ios::app | std::ios::out | std::ios::binary);
#endif
}

Object::~Object()
{
    for (std::map<int, Object*>::iterator sIt = mCommitted.mSubs.begin();
        sIt != mCommitted.mSubs.end();
        ++sIt)
    {
        delete sIt->second;
    }

    mOfile.close();
}

bool Object::WriteStep(const ObjectMgr* mgr, bool last)
{
    if (mOfile.is_open())
    {
        // save if time is less than orbitTime and we're saving this orbit,
        // given that we're skipping 10 orbitTimes for every one we save
//        double tstop = 0;
//        double tend = 0;
//        mgr->GetRunTimeData(tstop, tend);
//        double orbitTime = mgr->GetOrbitTime();
//        int norbits = (tstop - tend) / orbitTime;
//        norbits %= maxskip;
        int printSize = 10000;

        // need to save when (norbits == 0)
        // need to save and write when we reach (norbits == 1), then clear data
//        if ((norbits == 0) || ((norbits == 1) && !mOutput.empty()) || last)
        {
            mOutput.push_back(OutputLoc(mCommitted.mLoc[0], mCommitted.mLoc[1], mCommitted.mLoc[2]));

            // only write if we have a lot of data
            if ((mOutput.size() > printSize) || last)
//            if (last || (norbits == 1))
            {
                double step = mgr->GetStepTime();
                int buflen = 3 * mOutput.size() * sizeof(float);
                char* buf = new char[buflen];
                int len = 0;

                for (std::vector<OutputLoc>::iterator oIt = mOutput.begin();
                    oIt != mOutput.end();
                    ++oIt)
                {
                    float* fbuf = (float*)(buf + len);
                    fbuf[0] = oIt->mX;
                    fbuf[1] = oIt->mY;
                    fbuf[2] = oIt->mZ;
                    len += 3 * sizeof(float);
                }

                if (last)
                {
                    printf("End of run\n");
                }

                // write the data to file and clear it from the container
                mOfile.write(buf, len);
                delete[] buf;
                buf = 0;
                mOutput.clear();
            }
        }

        return true;
    }

    return false;
}

void Object::AddDeltaV(double dvx, double dvy, double dvz, Object* obj)
{
    DeltaV dv(dvx, dvy, dvz, obj);
    mCurrent.mDeltaV[obj->GetIndex()] = dv;
}

void Object::AppendDeltaV(double dvx, double dvy, double dvz, Object* obj)
{
    // create or replace
    std::map<int, DeltaV>::iterator dIt = mCurrent.mDeltaV.find(obj->GetIndex());

    if (dIt != mCurrent.mDeltaV.end())
    {
        // reset values to new entries
        DeltaV& dv = dIt->second;
        dv.mDvx = dvx;
        dv.mDvy = dvy;
        dv.mDvz = dvz;
    }
    else
    {
        AddDeltaV(dvx, dvy, dvz, obj);
    }
}

void Object::ProcessDeltaVs(double timestep)
{
    if (mCurrent.mDeltaV.empty())
    {
        // already handled by collision code
        return;
    }

    // initialize to committed values
    mCurrent.mLoc[0] = mCommitted.mLoc[0];
    mCurrent.mLoc[1] = mCommitted.mLoc[1];
    mCurrent.mLoc[2] = mCommitted.mLoc[2];
    mCurrent.mVx = mCommitted.mVx;
    mCurrent.mVy = mCommitted.mVy;
    mCurrent.mVz = mCommitted.mVz;

    // process triplets
    for (std::map<int, DeltaV>::iterator dIt = mCurrent.mDeltaV.begin();
        dIt != mCurrent.mDeltaV.end();
        ++dIt)
    {
        mCurrent.mVx += dIt->second.mDvx;
        mCurrent.mVy += dIt->second.mDvy;
        mCurrent.mVz += dIt->second.mDvz;
    }

    mCurrent.mLoc[0] += timestep * mCurrent.mVx;
    mCurrent.mLoc[1] += timestep * mCurrent.mVy;
    mCurrent.mLoc[2] += timestep * mCurrent.mVz;

    for (std::map<int, DeltaV>::iterator dIt = mCurrent.mDeltaV.begin();
        dIt != mCurrent.mDeltaV.end();
        ++dIt)
    {
        mCurrent.mVx += dIt->second.mDvx;
        mCurrent.mVy += dIt->second.mDvy;
        mCurrent.mVz += dIt->second.mDvz;
    }

    mCurrent.mDeltaV.clear();
}

void Object::ProcessFractionalDeltaVs(std::set<Object*>& multiGp, double timestep, double factor)
{
    std::set<int> multiIds;

    for (std::set<Object*>::iterator mIt = multiGp.begin();
        mIt != multiGp.end();
        ++mIt)
    {
        Object* obj = *mIt;
        multiIds.insert(obj->GetIndex());
    }

    // process triplets
    // count data from Objects NOT in multigp at reduced rate, 1/factor
    // take full value for objects IN multigp
    for (std::map<int, DeltaV>::iterator dIt = mCurrent.mDeltaV.begin();
        dIt != mCurrent.mDeltaV.end();
        ++dIt)
    {
        std::set<int>::iterator mIt = multiIds.find(dIt->first);

        // if one of the multiGp, take entirely
        if (mIt != multiIds.end())
        {
            mCurrent.mVx += dIt->second.mDvx;
            mCurrent.mVy += dIt->second.mDvy;
            mCurrent.mVz += dIt->second.mDvz;
        }
        else
        {
            mCurrent.mVx += dIt->second.mDvx / factor;
            mCurrent.mVy += dIt->second.mDvy / factor;
            mCurrent.mVz += dIt->second.mDvz / factor;
        }
    }

    mCurrent.mLoc[0] += timestep * mCurrent.mVx;
    mCurrent.mLoc[1] += timestep * mCurrent.mVy;
    mCurrent.mLoc[2] += timestep * mCurrent.mVz;

    for (std::map<int, DeltaV>::iterator dIt = mCurrent.mDeltaV.begin();
        dIt != mCurrent.mDeltaV.end();
        ++dIt)
    {
        std::set<int>::iterator mIt = multiIds.find(dIt->first);

        // if one of the multiGp, take entirely
        if (mIt != multiIds.end())
        {
            mCurrent.mVx += dIt->second.mDvx;
            mCurrent.mVy += dIt->second.mDvy;
            mCurrent.mVz += dIt->second.mDvz;
        }
        else
        {
            mCurrent.mVx += dIt->second.mDvx / factor;
            mCurrent.mVy += dIt->second.mDvy / factor;
            mCurrent.mVz += dIt->second.mDvz / factor;
        }
    }
}

bool Object::CalcLocFlags(double xmin, double xmax,
                          double ymin, double ymax,
                          double zmin, double zmax)
{
    // object xyz valuea
    // use these and offset to handle any cases where objects
    // are side by side but separated by a bit boundary
    double x = mCommitted.mLoc[0];
    double y = mCommitted.mLoc[1];
    double z = mCommitted.mLoc[2];

    // divide space into 32 chunks each for coarse, mid and fine
    double coarse = (xmax - xmin);

    if ((ymax - ymin) > coarse)
    {
        coarse = (ymax - ymin);
    }

    if ((zmax - zmin) > coarse)
    {
        coarse = (zmax - zmin);
    }

    // divide coarse into 32 chunks for mid
    double mid = coarse / 32;

    // divide mid into 32 chunks for fine
    double fine = mid / 32;

    mLocFlagsX.mCoarse = 1 << (int)(((x - xmin) / coarse) * 32.0);
    mLocFlagsY.mCoarse = 1 << (int)(((y - ymin) / coarse) * 32.0);
    mLocFlagsZ.mCoarse = 1 << (int)(((z - zmin) / coarse) * 32.0);

    mLocFlagsX.mMid = 1 << (int)(((x - (mLocFlagsX.mCoarse * mid)) / mid) * 32.0);
    mLocFlagsY.mMid = 1 << (int)(((y - (mLocFlagsY.mCoarse * mid)) / mid) * 32.0);
    mLocFlagsZ.mMid = 1 << (int)(((z - (mLocFlagsZ.mCoarse * mid)) / mid) * 32.0);

    mLocFlagsX.mFine = 1 << (int)(((x - (mLocFlagsX.mCoarse * mid) - (mLocFlagsX.mMid * fine)) / fine) * 32.0);
    mLocFlagsY.mFine = 1 << (int)(((y - (mLocFlagsY.mCoarse * mid) - (mLocFlagsY.mMid * fine)) / fine) * 32.0);
    mLocFlagsZ.mFine = 1 << (int)(((z - (mLocFlagsZ.mCoarse * mid) - (mLocFlagsZ.mMid * fine)) / fine) * 32.0);

    // object offset xyz values
    double xo = x + (mid / 2) + (fine / 2);
    double yo = y + (mid / 2) + (fine / 2);
    double zo = z + (mid / 2) + (fine / 2);

    mLocFlags2X.mCoarse = 1 << (int)(((xo - xmin) / coarse) * 32.0);
    mLocFlags2Y.mCoarse = 1 << (int)(((yo - ymin) / coarse) * 32.0);
    mLocFlags2Z.mCoarse = 1 << (int)(((zo - zmin) / coarse) * 32.0);

    mLocFlags2X.mMid = 1 << (int)(((xo - (mLocFlags2X.mCoarse * mid)) / mid) * 32.0);
    mLocFlags2Y.mMid = 1 << (int)(((yo - (mLocFlags2Y.mCoarse * mid)) / mid) * 32.0);
    mLocFlags2Z.mMid = 1 << (int)(((zo - (mLocFlags2Z.mCoarse * mid)) / mid) * 32.0);

    mLocFlags2X.mFine = 1 << (int)(((xo - (mLocFlags2X.mCoarse * mid) - (mLocFlags2X.mMid * fine)) / fine) * 32.0);
    mLocFlags2Y.mFine = 1 << (int)(((yo - (mLocFlags2Y.mCoarse * mid) - (mLocFlags2Y.mMid * fine)) / fine) * 32.0);
    mLocFlags2Z.mFine = 1 << (int)(((zo - (mLocFlags2Z.mCoarse * mid) - (mLocFlags2Z.mMid * fine)) / fine) * 32.0);

    //printf("Setting locflags, object %d :  coarse flags %x, %x, %x, and half over: %x, %x, %x\n",
    //    GetIndex(),
    //    GetLocFlagsX().mCoarse,
    //    GetLocFlagsY().mCoarse,
    //    GetLocFlagsZ().mCoarse,
    //    GetLocFlags2X().mCoarse,
    //    GetLocFlags2Y().mCoarse,
    //    GetLocFlags2Z().mCoarse);

    return true;
}

//bool Object::IntersectsInTimestep(Object* obj, 
//                                  double timestep, 
//                                  double& proximityTime0,
//                                  double& proximityTime1,
//                                  double& sqdist,
//                                  double closeFactor,
//                                  int currentTimestep)
//{
//    // if there's a locflags overlap, check distance between travel vectors
//    if (((obj->mLocFlagsX.mCoarse & mLocFlagsX.mCoarse) ||
//        (obj->mLocFlags2X.mCoarse & mLocFlags2X.mCoarse)) &&
//        ((obj->mLocFlagsY.mCoarse & mLocFlagsY.mCoarse) ||
//        (obj->mLocFlags2Y.mCoarse & mLocFlags2Y.mCoarse)) &&
//        ((obj->mLocFlagsZ.mCoarse & mLocFlagsZ.mCoarse) ||
//        (obj->mLocFlags2Z.mCoarse & mLocFlags2Z.mCoarse)))
//    {
////        printf("Flag overlap between ID %d and ID %d, num steps = %d\n", mIndex, obj->GetIndex(), currentTimestep);
//        // check whether there's any proximity between the objects during this timestep
//        double x3 = 0;
//        double y3 = 0;
//        double z3 = 0;
//        double vx2 = 0;
//        double vy2 = 0;
//        double vz2 = 0;
//
//        obj->GetLocation(x3, y3, z3);
//        obj->GetVelocity(vx2, vy2, vz2);
//
//        double x4 = x3 + (timestep * vx2);
//        double y4 = y3 + (timestep * vy2);
//        double z4 = z3 + (timestep * vz2);
//
//        double x1 = mCommitted.mLoc[0];
//        double y1 = mCommitted.mLoc[1];
//        double z1 = mCommitted.mLoc[2];
//
//        double x2 = x1 + (timestep * mCommitted.mVx);
//        double y2 = y1 + (timestep * mCommitted.mVy);
//        double z2 = z1 + (timestep * mCommitted.mVz);
//
//        double rad2 = obj->GetRadius();
//
//        // if we want to know precisely whether objects are colliding, use exact
//        double bothrad = rad2 + mRadius;
//
//        // The following algorithm posted online without restriction by Paul Bourke, dated April 1998
//        //
//        // a line will be defined by two points lying on it, 
//        // a point on line "a" defined by points P1 and P2 has an equation.
//        //    Pa = P1 + mua(P2 - P1)
//        // similarly a point on a second line "b" defined by points P4 and P4 will be written as
//        //    Pb = P3 + mub(P4 - P3)
//        //
//        // The values of mua and mub range from negative to positive infinity.The line segments between P1 P2 and P3 P4 have their corresponding mu between 0 and 1.
//        //
//        // The shortest line segment between two 3D lines will be perpendicular to the two lines.
//        // This allows us to write two equations for the dot product as
//        //    (Pa - Pb) dot(P2 - P1) = 0
//        //    (Pa - Pb) dot(P4 - P3) = 0
//        // where Pa and Pb are the intersections of normals on the lines 
//        //
//        // Expanding these given the equation of the lines
//        //    (P1 - P3 + mua(P2 - P1) - mub(P4 - P3)) dot(P2 - P1) = 0
//        //    (P1 - P3 + mua(P2 - P1) - mub(P4 - P3)) dot(P4 - P3) = 0
//        //
//        //    Expanding these in terms of the coordinates(x, y, z) is a nightmare but the result is as follows
//        //
//        //    d1321 + mua d2121 - mub d4321 = 0
//        //    d1343 + mua d4321 - mub d4343 = 0
//        //
//        // where
//        //    dmnop = (xm - xn)(xo - xp) + (ym - yn)(yo - yp) + (zm - zn)(zo - zp)
//        //    Note that dmnop = dopmn
//        //
//        //    Finally, solving for mua gives
//        //
//        //    mua = (d1343 d4321 - d1321 d4343) / (d2121 d4343 - d4321 d4321)
//        //    and back - substituting gives mub
//        //
//        //    mub = (d1343 + mua d4321) / d4343
//        // 
//        double xd13 = (x1 - x3);
//        double xd43 = (x4 - x3);
//        double xd21 = (x2 - x1);
//        double yd13 = (y1 - y3);
//        double yd43 = (y4 - y3);
//        double yd21 = (y2 - y1);
//        double zd13 = (z1 - z3);
//        double zd43 = (z4 - z3);
//        double zd21 = (z2 - z1);
//        double d1343 = xd13 * xd43 + yd13 * yd43 + zd13 * zd43;
//        double d4321 = xd43 * xd21 + yd43 * yd21 + zd43 * zd21;
//        double d1321 = xd13 * xd21 + yd13 * yd21 + zd13 * zd21;
//        double d4343 = xd43 * xd43 + yd43 * yd43 + zd43 * zd43;
//        double d2121 = xd21 * xd21 + yd21 * yd21 + zd21 * zd21;
//        double mua = (d1343 * d4321 - d1321 * d4343) / (d2121 * d4343 - d4321 * d4321);
//        double mub = (d1343 + mua * d4321) / d4343;
//        double pax = x1 + mua * (x2 - x1);
//        double pay = y1 + mua * (y2 - y1);
//        double paz = z1 + mua * (z2 - z1);
//
//        // make sure point is not beyond line segment
//        if (((pax < x1) && (x1 < x2)) || ((pax > x1) && (x1 > x2)))
//        {
//            pax = x1;
//        }
//        else if (((pax > x2) && (x1 < x2)) || ((pax < x2) && (x1 > x2)))
//        {
//            pax = x2;
//        }
//
//        if (((pay < y1) && (y1 < y2)) || ((pay > y1) && (y1 > y2)))
//        {
//            pay = y1;
//        }
//        else if (((pay > y2) && (y1 < y2)) || ((pay < y2) && (y1 > y2)))
//        {
//            pay = y2;
//        }
//
//        if (((paz < z1) && (z1 < z2)) || ((paz > z1) && (z1 > z2)))
//        {
//            paz = z1;
//        }
//        else if (((paz > z2) && (z1 < z2)) || ((paz < z2) && (z1 > z2)))
//        {
//            paz = z2;
//        }
//
//        double pbx = x3 + mub * (x4 - x3);
//        double pby = y3 + mub * (y4 - y3);
//        double pbz = z3 + mub * (z4 - z3);
//
//        // make sure other point is not beyond line segment
//        if (((pbx < x3) && (x3 < x4)) || ((pbx > x3) && (x3 > x4)))
//        {
//            pbx = x3;
//        }
//        else if (((pbx > x4) && (x3 < x4)) || ((pbx < x4) && (x3 > x4)))
//        {
//            pbx = x4;
//        }
//
//        if (((pby < y3) && (y3 < y4)) || ((pby > y3) && (y3 > y4)))
//        {
//            pby = y3;
//        }
//        else if (((pby > y4) && (y3 < y4)) || ((pby < y4) && (y3 > y4)))
//        {
//            pby = y4;
//        }
//
//        if (((pbz < z3) && (z3 < z4)) || ((pbz > z3) && (z3 > z4)))
//        {
//            pbz = z3;
//        }
//        else if (((pbz > z4) && (z3 < z4)) || ((pbz < z4) && (z3 > z4)))
//        {
//            pbz = z4;
//        }
//
//        double xab = pax - pbx;
//        double yab = pay - pby;
//        double zab = paz - pbz;
//        double rsq = (xab * xab) + (yab * yab) + (zab * zab);
//        sqdist = rsq;
//
//        // check whether the distance between the traversal vectors is less than 
//        // the combined radii (or twice that if "exact" not specified)
//        if (rsq < (bothrad * bothrad * closeFactor * closeFactor))
//        {
//            // calculate the times of closest approach
//            proximityTime0 = abs((pax - x1) / mCommitted.mVx);
//            double pty = abs((pay - y1) / mCommitted.mVy);
//            double ptz = abs((paz - z1) / mCommitted.mVz);
//
//            if (pty > proximityTime0)
//            {
//                proximityTime0 = pty;
//            }
//
//            if (ptz > proximityTime0)
//            {
//                proximityTime0 = ptz;
//            }
//
//            proximityTime1 = abs((pbx - x3) / vx2);
//
//            double pty1 = abs((pby - y3) / vy2);
//            double ptz1 = abs((pbz - z3) / vz2);
//
//            if (pty1 > proximityTime1)
//            {
//                proximityTime1 = pty1;
//            }
//
//            if (ptz1 > proximityTime1)
//            {
//                proximityTime1 = ptz1;
//            }
//
//            return true;
//        }
//    }
//    else
//    {
////        printf("Skipping distance check for ID %d and ID %d\n", mIndex, obj->GetIndex());
//    }
//
//    return false;
//}

void Object::RemoveDeltaVandInitialize(Object* obj)
{
    // remove and make sure we;re initialized, as well
    mCurrent.mLoc[0] = mCommitted.mLoc[0];
    mCurrent.mLoc[1] = mCommitted.mLoc[1];
    mCurrent.mLoc[2] = mCommitted.mLoc[2];
    mCurrent.mVx = mCommitted.mVx;
    mCurrent.mVy = mCommitted.mVy;
    mCurrent.mVz = mCommitted.mVz;

    std::map<int, DeltaV>::iterator dIt = mCurrent.mDeltaV.find(obj->GetIndex());

    if (dIt != mCurrent.mDeltaV.end())
    {
        mCurrent.mDeltaV.erase(dIt);
    }
}

bool Object::IntersectsInTimestep(Object* obj,
                                  double& proximityTime,
                                  double& dist,
                                  double closeFactor,
                                  double currentTimestep)
{
    // if there's a locflags overlap, check distance between travel vectors
    if (((obj->mLocFlagsX.mCoarse & mLocFlagsX.mCoarse) ||
        (obj->mLocFlags2X.mCoarse & mLocFlags2X.mCoarse)) &&
        ((obj->mLocFlagsY.mCoarse & mLocFlagsY.mCoarse) ||
        (obj->mLocFlags2Y.mCoarse & mLocFlags2Y.mCoarse)) &&
            ((obj->mLocFlagsZ.mCoarse & mLocFlagsZ.mCoarse) ||
        (obj->mLocFlags2Z.mCoarse & mLocFlags2Z.mCoarse)))
    {
        // NOTE:  ignore contact at beginning of timestep

        // check whether there's any proximity between the objects during this timestep
        // D0^2 = dx^2 + dy^2 + dz^2
        // D1^2 = dx^2 + dy^2 + dz^2 + 2t(dx * dvx + dy * dvy + dz * dvz) + t^2(dvx^2 + dvy^2 + dvz^2)
        // to find the min for D1^2 during time t, differentiate and set to 0
        // D0^2 - D1^2 = -2t(dx * dvx + dy * dvy + dz * dvz) - t^2(dvx^2 + dvy^2 + dvz^2)
        // d/dt of D1^2 (with dv's as constants) = 2(dx * dvx + dy * dvy + dz * dvz) + 2t(dvx^2 + dvy^2 + dvz^2) = 0
        // => (min or max dist) t = -(dx * dvx + dy * dvy + dz * dvz) / (dvx^2 + dvy^2 + dvz^2)
        double x0 = 0;
        double y0 = 0;
        double z0 = 0;
        GetLocation(x0, y0, z0);
        double vx0 = 0;
        double vy0 = 0;
        double vz0 = 0;
        GetVelocity(vx0, vy0, vz0);

        double x1 = 0;
        double y1 = 0;
        double z1 = 0;
        obj->GetLocation(x1, y1, z1);
        double vx1 = 0;
        double vy1 = 0;
        double vz1 = 0;
        obj->GetVelocity(vx1, vy1, vz1);

        double dx = x0 - x1;
        double dy = y0 - y1;
        double dz = z0 - z1;

        double dvx = vx0 - vx1;
        double dvy = vy0 - vy1;
        double dvz = vz0 - vz1;

        double dxdv = dx * dvx;
        double dydv = dy * dvy;
        double dzdv = dz * dvz;

        double dvx2 = dvx * dvx;
        double dvy2 = dvy * dvy;
        double dvz2 = dvz * dvz;

        double rad2 = obj->GetRadius();
        double bothrad = rad2 + mRadius;

        double distS2 = dx * dx + dy * dy + dz * dz;
        double mt = -(dx * dvx + dy * dvy + dz * dvz) / (dvx * dvx + dvy * dvy + dvz * dvz);

        double distS = sqrt(distS2);
        // note: initial distance should never be less that surfaces touching
//        assert(distS >= bothrad);

        double distMT2 = distS2;
        double distMT = distS;
        double distE2 = distS2 + 2 * currentTimestep * (dxdv + dydv + dzdv) +
            currentTimestep * currentTimestep * (dvx2 + dvy2 + dvz2);
        double distE = sqrt(distE2);

        // if we're not really that close yet and collision time is a ways out
        // return false for now
        if ((distE > (closeFactor * bothrad)) && ((mt < 0) || (mt > currentTimestep)))
        {
            // not in range
            return false;
        }

        double tstep = 0;
        bool minS = true;

        // set up so distMT is the min dist and tstep is appropriate
        // for distS closest, set above
        if ((mt > 0) && (mt < currentTimestep))
        {
            distMT2 = distS2 + 2 * mt * (dxdv + dydv + dzdv) + mt * mt* (dvx2 + dvy2 + dvz2);
            distMT = sqrt(distMT2);
        }

        // if we have a minimum in this timestep, set tstep
        if ((distMT <= distE) && (distMT < distS) && (mt > 0) && (mt < currentTimestep))
        {
            // use mt to calculate proximityTime
            tstep = mt;
            minS = false;
        }
        else if ((distE < distMT) && (distE < distS))
        {
            distMT = distE;
            // use mt to calculate proximityTime
            tstep = currentTimestep;
            minS = false;
        }

        // check beginning and mt
        if (!minS && (distMT <= bothrad))
        {
            // collision in this timestep
            // approach is linear (as a first approximation), so get time from proportion
            // note that we can't have a collision at beginning, so distMT should never be distS
            proximityTime = tstep * (distS - bothrad) / (distS - distMT);
            dist = bothrad;
            return true;
        }
        else if (!minS && (distMT < (closeFactor * bothrad)))
        {
            proximityTime = tstep;
            dist = distMT;
            return true;
        }
    }

    return false;
}

double Object::CheckTimestepCriteria(void) const
{
    double tsfactor = 1;

    // criteria for timestep adjustment:
    // 1 - rate of momentum change
    //    calculate |dv|^2 / |v|^2 (unless v is close to 0), use as a criterion
    // 2 - rate of force angle change as a proportion of velocity
    //    calculate |v|^2 - (dv dot v)^2 / |v|^2 = dv(normal)^2
    //    use dv(normal)^2 / |v|^2 as a criterion

    double magV2 = mCommitted.MagV2();
    double dvx = 0;
    double dvy = 0;
    double dvz = 0;
    double magDV2 = mCurrent.MagDV2(dvx, dvy, dvz);

    if (abs(magV2) > kMinVforDVCheck)
    {
        double dmom = magDV2 / magV2;

        // if dmom is bigger than the square of our tolerance, 
        // set tsfactor large enough to get the value in bounds
        if (dmom > (kMaxDVChange * kMaxDVChange))
        {
            // use sqrt to determine tsfactor
            tsfactor = (int)(sqrt(dmom) / kMaxDVChange) + 1;
        }
    }

    double dot = dvx * mCommitted.mVx + dvy * mCommitted.mVy + dvz * mCommitted.mVz;
    double dnorm = abs(1 - dot * dot / (magV2 * magV2));

    if (dnorm > (kMaxDVNormalChange * kMaxDVNormalChange))
    {
        double factor = (int)(sqrt(dnorm) / kMaxDVNormalChange) + 1;

        if (factor > tsfactor)
        {
            tsfactor = factor;
        }
    }

    return tsfactor;
}

double TimeStepData::MagV2(void) const
{
    return mVx * mVx + mVy * mVy + mVz * mVz;
}

double TimeStepData::MagDV2(double& dvx, double& dvy, double& dvz) const
{
    for (std::map<int, DeltaV>::const_iterator dIt = mDeltaV.begin();
        dIt != mDeltaV.end();
        ++dIt)
    {
        dvx += dIt->second.mDvx;
        dvy += dIt->second.mDvy;
        dvz += dIt->second.mDvz;
    }

    return dvx * dvx + dvy * dvy + dvz * dvz;
}
