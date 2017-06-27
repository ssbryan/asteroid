#include "ObjectMgr.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <set>
#include <ctime>
#include <assert.h>

// Author Stephen Bryan
// June 22, 2017

void FindCollisionGroups(std::vector<CollisionData>& collisions, std::vector<std::vector<CollisionData> >& multiCollisions);

ObjectMgr::ObjectMgr(const std::string& fname)
    : mTEnd(1e6)
    , mTStep(1)
    , mTStepOriginal(1)
    , mTStop(1e6)
    , mPrint(0)
    , mMaxOrbitTime(1)
    , mXmax(0)
    , mXmin(0)
    , mYmax(0)
    , mYmin(0)
    , mZmax(0)
    , mZmin(0)
{
    mOkay = true;
    Initialize(fname);
}

bool ObjectMgr::Initialize(const std::string& fname)
{
    time(&mStartTime);
    const char* fstr = fname.c_str();
    std::string logExt(".log");

    if (*(fstr + strlen(fstr) - 1) == '.')
    {
        logExt = "log";
    }

    std::string logname = fname + logExt;
#ifndef PC
        mLogFile.open(logname.c_str(), std::ios_base::out);
#else
        mLogFile.open(logname.c_str(), std::ios::out);
#endif

    std::ifstream datafile;
    int numobjects = 0;
    double maxradius = 1;

    // open the named file, read in initialization data and create objects
    // format of file (all doubles):  x y z vx vy vz mass\n
#ifndef PC
    datafile.open(fname.c_str(), std::ios_base::in | std::ios_base::binary);
#else
    datafile.open(fname, std::ios::in | std::ios::binary);
#endif

    char msgbuf4[80];
    sprintf(msgbuf4, "Initialize:  Reading data file %s\n", fname.c_str());
    LogSimpleMessage(msgbuf4);

    if (datafile.is_open())
    {
        bool okay = true;
        // format for first line: Options: opt1name-number, opt2name-number,...\n
        // recognized options with defaults are:
        // momtol=1e-2
        // dvnormtol=1e-2
        // runs=1e9
        // tstep=1
        // minvforcheck=1.0
        // startdatasave=0
        // nthdatasave=1

        char options[200];
        datafile.getline(options, 199);

        if (strnicmp(options, "Options:", 8))
        {
            // no Options line
            mOkay = false;
            std::string msg("Initialize:  Initialization failed - no Options line in file ");
            msg += fname + "\n";
            LogSimpleMessage(msg.c_str());
            datafile.close();
            return false;
        }

        // record datafile info
        LogSimpleMessage(options);
        LogSimpleMessage("\n");
        char name[40];
        char val[40];
        char* opt = name;

        for (char* c = options + 8; ; ++c)
        {
            if ((*c == ' ') || (*c == '\t'))
            {
                // skip
                continue;
            }

            if ((*c == ',') || (*c == 0))
            {
                // finished with option:  set value
                *opt = 0;
                char* endptr = 0;

                if (!stricmp(name, "momtol"))
                {
                    kMaxDVChange = strtod(val, &endptr);
                }
                else if (!stricmp(name, "dvnormtol"))
                {
                    kMaxDVNormalChange = strtod(val, &endptr);
                }
                else if (!stricmp(name, "runs"))
                {
                    mTEnd = strtod(val, &endptr);
                    mTStop = mTEnd;
                }
                else if (!stricmp(name, "tstep"))
                {
                    mTStep = strtod(val, &endptr);
                    mTStepOriginal = mTStep;
                }
                else if (!stricmp(name, "minvforcheck"))
                {
                    kMinVforDVCheck = strtod(val, &endptr);
                }
                else if (!stricmp(name, "startdatasave"))
                {
                    kStartDataSave = strtod(val, &endptr);
                }
                else if (!stricmp(name, "nthdatasave"))
                {
                    // set value to 1 or greater
                    kNthDataSave = strtoul(val, &endptr, 10);

                    if (kNthDataSave == 0)
                    {
                        kNthDataSave = 1;
                    }
                }

                // and clear bufs
                name[0] = 0;
                val[0] = 0;
                opt = name;

                if (*c == 0)
                {
                    // end of line
                    break;
                }
            }
            else if (*c == '=')
            {
                *opt = 0;
                opt = val;
            }
            else
            {
                *opt++ = *c;
            }
        }

        while (okay)
        {
            if (datafile.eof())
            {
                double totmass = 0;
                double comx = 0;
                double comy = 0;
                double comz = 0;
                double pvx = 0;
                double pvy = 0;
                double pvz = 0;

                // get center of mass, total momentum
                for (std::vector<Object*>::iterator oIt = mObjs.begin();
                    oIt != mObjs.end();
                    ++oIt)
                {
                    Object* obj = *oIt;
                    double x = 0;
                    double y = 0;
                    double z = 0;
                    obj->GetLocation(x, y, z);
                    double vx = 0;
                    double vy = 0;
                    double vz = 0;
                    obj->GetVelocity(vx, vy, vz);
                    double mass = obj->GetMass();
                    totmass += mass;
                    comx += x * mass;
                    comy += y * mass;
                    comz += z * mass;
                    pvx += vx * mass;
                    pvy += vy * mass;
                    pvz += vz * mass;
                }

                char msgbuf[80];
                sprintf(msgbuf, "Initialize:  Initialization completed.  %d particles read in\n", numobjects);
                LogSimpleMessage(msgbuf);
                comx /= totmass;
                comy /= totmass;
                comz /= totmass;
                double angmx = 0;
                double angmy = 0;
                double angmz = 0;

                // determine angular momentum
                for (std::vector<Object*>::iterator oIt = mObjs.begin();
                    oIt != mObjs.end();
                    ++oIt)
                {
                    Object* obj = *oIt;
                    double x = 0;
                    double y = 0;
                    double z = 0;
                    obj->GetLocation(x, y, z);
                    double vx = 0;
                    double vy = 0;
                    double vz = 0;
                    obj->GetVelocity(vx, vy, vz);
                    double mass = obj->GetMass();
                    angmx += (comx - x) * vx * mass;
                    angmy += (comy - y) * vy * mass;
                    angmz += (comz - z) * vz * mass;

                    double r = sqrt((comx - x) * (comx - x) +
                        (comy - y) * (comy - y) +
                        (comz - z) * (comz - z));
                    double gmass = obj->GetGMass();
                    double cmass = totmass - gmass;
                    // here's the ideal orbital velocity at this location and this mass
                    // ***assuming*** a large central object is most of the mass in the system
                    double orbitalVel = (r < 1e-6) ? 0 : sqrt(cmass * G / r);

                    if (orbitalVel > 0)
                    {
                        double orbitTime = 2 * r * 3.14159265 / orbitalVel;

                        if (orbitTime > mMaxOrbitTime)
                        {
                            mMaxOrbitTime = orbitTime;
                        }
                    }

                    // now get the actual perpendicular (to the vector to COM) and parallel velocities
                    double vtot2 = vx * vx + vy * vy + vz * vz;
                    double rz2 = sqrt(r * r - (comz - z) * (comz - z));
                    double vz2 = (vtot2 < 1e-6) ? 0 : sqrt(1 - (vz * vz) / vtot2);
                    bool ctr = rz2 < 1e-6;
                    double vpe = ctr ? 0 : ((x - comx) * vy + (comy - y) * vx) / rz2 * vz2;
                    double vpa = ctr ? 0 : ((comy - y) * vy + (comx - x) * vx) / rz2 * vz2;

                    char msgbuf1[120];
                    sprintf(msgbuf1, "Initialize:  Object %d, ideal orbital vel: \t%g; perpendicular: \t%g, parallel: \t%g\n",
                        obj->GetIndex(), orbitalVel, vpe, vpa);
                    LogSimpleMessage(msgbuf1);
                }

                char msgbuf2[120];
                sprintf(msgbuf2, "Initialize:  Center of mass: \t%g, %g, %g\nMomentum: \t\t%g, %g, %g\nAngular momentum: \t%g, %g, %g\nMaxOrbitTime: %g\n",
                    comx, comy, comz, pvx, pvy, pvz, angmx, angmy, angmz, mMaxOrbitTime);
                LogSimpleMessage(msgbuf2);
                okay = false;
                continue;
            }

            double x = 0;
            double y = 0;
            double z = 0;
            double vx = 0;
            double vy = 0;
            double vz = 0;
            double mass = 0;
            datafile >> x >> y >> z >> vx >> vy >> vz >> mass;

            char msgbuf0[120];
            sprintf(msgbuf0, "%g %g %g %g %g %g %g\n",
                x, y, z, vx, vy, vz, mass);
            LogSimpleMessage(msgbuf0);

            if (datafile.eof())
            {
                continue;
            }
            else if (datafile.fail())
            {
                okay = false;
                mOkay = false;
                LogSimpleMessage("Initialize:  Initialization failed.\n");
            }
            else
            {
                Object* obj = new Object(x, y, z, vx, vy, vz, mass, numobjects);

                if (x > mXmax)
                {
                    mXmax = x;
                }

                if (y > mYmax)
                {
                    mYmax = y;
                }

                if (z > mZmax)
                {
                    mZmax = z;
                }

                if (x < mXmin)
                {
                    mXmin = x;
                }

                if (y < mYmin)
                {
                    mYmin = y;
                }

                if (z < mZmin)
                {
                    mZmin = z;
                }

                double rad = obj->GetRadius();

                if (rad > maxradius)
                {
                    maxradius = rad;
                }

                mObjs.push_back(obj);
                numobjects++;
            }
        }

        datafile.close();
    }

    mMaxOrbitTime *= 2;

    // increase the extent in all directions
    double diff = mXmax - mXmin;

    if ((mYmax - mYmin) > diff)
    {
        diff = mYmax - mYmin;
    }

    if ((mZmax - mZmin) > diff)
    {
        diff = mZmax - mZmin;
    }

    // set minimum coarseness so objects are in the same coarse box when close
    if (diff < (4 * 32 * maxradius))
    {
        diff = 4 * 32 * maxradius;
    }

    mXmin = (mXmax + mXmin - diff) / 2;
    mXmax = mXmin + diff;
    mYmin = (mYmax + mYmin - diff) / 2;
    mYmax = mYmin + diff;
    mZmin = (mZmax + mZmin - diff) / 2;
    mZmax = mZmin + diff;

    char msgbuf3[120];
    sprintf(msgbuf3, "Initialize: Space extents: %g : %g, %g : %g, %g : %g\n",
        mXmin, mXmax, mYmin, mYmax, mZmin, mZmax);
    LogSimpleMessage(msgbuf3);

    //// set initial LocFlags
    //for (std::vector<Object*>::iterator oIt = mObjs.begin();
    //    oIt != mObjs.end();
    //    ++oIt)
    //{
    //    Object* obj = *oIt;
    //    obj->CalcLocFlags(diff, mXmin, mYmin, mZmin);
    //    printf("ObjectMgr::Initialize: Setting locflags, object %d :  coarse flags %x, %x, %x, and half over: %x, %x, %x\n",
    //        obj->GetIndex(),
    //        obj->GetLocFlagsX().mCoarse,
    //        obj->GetLocFlagsY().mCoarse,
    //        obj->GetLocFlagsZ().mCoarse,
    //        obj->GetLocFlags2X().mCoarse,
    //        obj->GetLocFlags2Y().mCoarse,
    //        obj->GetLocFlags2Z().mCoarse);
    //}

    return numobjects > 1;
}

bool ObjectMgr::Run(void)
{
    if (!mOkay)
    {
        return false;
    }

    // each pass, all Object flags are one way or the other
    // don't make a separate pass resetting them,
    // but just change the sense of what's "done"
    bool done = false;

    while (CalcStep(done))
    {
        done = !done;
    }

    // print/log the elapsed time
    int hrs = 0;
    int mins = 0;
    int secs = 0;
    time_t endTime = mStartTime;
    time(&endTime);
    double elapsed = difftime(endTime, mStartTime);

    while (elapsed > 3600)
    {
        elapsed -= 3600;
        hrs++;
    }

    while (elapsed > 60)
    {
        elapsed -= 60;
        mins++;
    }

    char msgbuf[80];
    sprintf(msgbuf, "\nRun:  Elapsed time: %d hrs, %d mins, %g secs",
        hrs, mins, elapsed);
    LogSimpleMessage(msgbuf);

    mLogFile.close();

    // pause
//    getchar();
    return true;
}


bool ObjectMgr::CalcStep(bool done)
{
    // running total (which may roll over) of steps, for saving data
    static unsigned int stepnum = 0;
    // count dots printed, use to output fractions occasionally
    static unsigned int dotnum = 0;

    // at the first step, run through once and caculate deltaVs
    // then get appropriate tsFactor, reset timestep, reset deltaVs, and begin for real
    if (stepnum == 0)
    {
        double tsFactor = 1;

        // first, calculate deltaVs fpr all objects
        for (std::vector<Object*>::iterator oIt = mObjs.begin();
            oIt != mObjs.end();
            ++oIt)
        {
            Object* obj = *oIt;

            for (std::vector<Object*>::iterator oIt2 = mObjs.begin();
                oIt2 != mObjs.end();
                ++oIt2)
            {
                Object* obj2 = *oIt2;

                if (obj == obj2)
                {
                    continue;
                }

                if (obj->IsDone(obj2, done))
                {
                    continue;
                }

                std::vector<CollisionData> multicolls;
                CollisionData cd(obj, obj2);
                multicolls.push_back(cd);

                CalcStepForMultiObjects(multicolls, mTStep);
                obj->SetDone(obj2, done);
                obj2->SetDone(obj, done);
            }
        }

        // then check constraints to get adjusted timestep
        for (std::vector<Object*>::iterator oIt = mObjs.begin();
            oIt != mObjs.end();
            ++oIt)
        {
            Object* obj = *oIt;
            double tsFactorTmp = obj->CheckTimestepCriteria();

            if (tsFactorTmp > tsFactor)
            {
                tsFactor = tsFactorTmp;
            }
        }

        // then reset deltaVs and done flags
        for (std::vector<Object*>::iterator oIt = mObjs.begin();
            oIt != mObjs.end();
            ++oIt)
        {
            Object* obj = *oIt;
            obj->ClearCurrentData();

            for (std::vector<Object*>::iterator oIt2 = mObjs.begin();
                oIt2 != mObjs.end();
                ++oIt2)
            {
                Object* obj2 = *oIt2;

                if (obj == obj2)
                {
                    continue;
                }

                if (obj->IsDone(obj2, !done))
                {
                    continue;
                }

                obj->SetDone(obj2, !done);
                obj2->SetDone(obj, !done);
            }
        }

        if (tsFactor > 1)
        {
            mTStep = mTStep / tsFactor;

            char msgbuf1[120];
            sprintf(msgbuf1, "\nCalcStep:  Initial calulation: TimeStep reset from %g to %g\n", mTStepOriginal, mTStep);
            LogSimpleMessage(msgbuf1);
        }
    }

    // record collisions
    std::vector<CollisionData> collisions;
    std::set<int> collidedObjects;

    // process Objects to get acceleration
    for (std::vector<Object*>::iterator oIt = mObjs.begin();
        oIt != mObjs.end();
        ++oIt)
    {
        Object* obj = *oIt;

        for (std::vector<Object*>::iterator oIt2 = mObjs.begin();
            oIt2 != mObjs.end();
            ++oIt2)
        {
            Object* obj2 = *oIt2;

            if (obj == obj2)
            {
                continue;
            }

            if (obj->IsDone(obj2, done))
            {
                continue;
            }

            // check for whether traversed volumes collide
            double tc = 0;
            double dist = 0;
            bool isClose = obj->IntersectsInTimestep(obj2, tc, dist, 2, mTStep);

            // don't look for exact collisions, just close (2 * combined radii)
            if (isClose)
            {
                // close encounter detected
                // save index pair so we can process later
                CollisionData cd(obj, obj2);
                cd.mDist = dist;
                cd.mTimeToCollision = tc;
                collisions.push_back(cd);
                collidedObjects.insert(obj->GetIndex());
                collidedObjects.insert(obj2->GetIndex());
            }

            // process anyway so we have a complete set of data in deltaVs and mark what's done consistently
            std::vector<CollisionData> multicolls;
            // NOTE: only Object* are valid in CollisionData below, since isClose is false
            CollisionData cd(obj, obj2);
            multicolls.push_back(cd);

            CalcStepForMultiObjects(multicolls, mTStep);
            obj->SetDone(obj2, done);
            obj2->SetDone(obj, done);
        }
    }

    // if any objects have come close, process in smaller timesteps
    if (!collisions.empty())
    {
        ProcessCollisions(collisions);
    }

    // tell objects the accelerations are all accumulated and it's time to process them
    for (std::vector<Object*>::iterator oIt = mObjs.begin();
        oIt != mObjs.end();
        ++oIt)
    {
        Object* obj = *oIt;
        bool inColliders = !collidedObjects.empty();

        if (inColliders)
        {
            std::set<int>::iterator cIt = collidedObjects.find(obj->GetIndex());
            inColliders = (cIt != collidedObjects.end());
        }

        if (!inColliders)
        {
            // process deltaVs and commit only if not involved in a collision
            obj->ProcessDeltaVs(mTStep);
            obj->CommitTimeStepData();
            obj->ClearCurrentData();
        }

        obj->CalcLocFlags(mXmax - mXmin, mXmin, mYmin, mZmin);

        // start writing data after optional kStartDataSave, defaulting to 0
        if ((mTStop - mTEnd) > kStartDataSave)
        {
            if ((kNthDataSave > 1) && (stepnum % kNthDataSave == 0))
            {
                obj->WriteStep(this, (mTEnd - mTStep) <= 0);
            }
        }
    }

    stepnum++;
    mTEnd -= mTStep;
    mPrint += mTStep;

    if (mPrint >= mTStop / 80)
    {
        dotnum++;

        if ((dotnum % 8) == 0)
        {
            // can't use LogSimpleMessage due to % in msg
            char msgbuf[120];
            sprintf(msgbuf, "  %d%%  ", (dotnum * 10) / 8);
            printf("  %d%%  ", (dotnum * 10) / 8);

            if (mLogFile.is_open())
            {
                mLogFile << msgbuf;
            }
        }

        printf(".");
        mPrint = 0;
    }

    return mTEnd > 0;
}

void ObjectMgr::ProcessCollisions(std::vector<CollisionData>& collisions)
{
    // NOTE:  need to process everything in the list - the Objects involved haven't been adjusted in this timestep

    // vector of groups of collision data that have potentially multiple objects 
    // all getting close to or colliding with one another
    std::vector<std::vector<CollisionData> > multiCollisions;
    FindCollisionGroups(collisions, multiCollisions);

    // now process the groups

    // for each group in the list
    // determine whether collision or close approach exists
    //   for close approach, 
    //      remove appropriate deltaVs
    //      calculate the new timestep
    //      apply forces for this group of objects to get new deltaVs
    //      check distance and make sure it's still reasonable
    //      after all mini timsteps, commit changes
    //   for collision,
    //      remove appropriate deltaVs
    //      calculate the new timestep to collision(s)
    //      move forward until collision is within one (mini) timestep
    //      calculate the change in momentum, velocity vector, location for 2 timesteps
    //        (one leading up to the collision, one to see where the objects are after collision)
    //      calculate the rest of the mini timesteps to the end of the main one
    //      commit changes

    // if distance is close but not changing much, 
    // just make sure we're not taking too coarse a timestep - may be okay as is
    for (std::vector<std::vector<CollisionData> >::iterator mIt = multiCollisions.begin();
        mIt != multiCollisions.end();
        ++mIt)
    {
        std::set<Object*> multiGp;

        // get min criterion (rollback, reduced timestep or pass)
        // for entire group
        bool rollback = false;
        double tsFactor = 1;

        for (std::vector<CollisionData>::iterator cIt = collisions.begin();
            cIt != collisions.end();
            ++cIt)
        {
            Object* obj0 = cIt->mObj0;
            Object* obj1 = cIt->mObj1;
            // we want to process all objects in this collision set
            // so accumulate unique elements here
            multiGp.insert(obj0);
            multiGp.insert(obj1);
            double dist = cIt->mDist;
            assert(dist > 0);

            // combined radii
            double rad2 = obj0->GetRadius() + obj1->GetRadius();
            double tol = 1.01;

            // if very close to actual collision, examine in detail
            if (dist < (tol * rad2))
            {
                char msgbuf[120];
                sprintf(msgbuf, "\nProcessCollisions:  Distance between ID %d and ID %d is %g (radius both = %g), time %g",
                    obj0->GetIndex(), obj1->GetIndex(), dist - rad2, rad2, mTStop - mTEnd);
                LogSimpleMessage(msgbuf);

                // distance is less than combined radii between the objects
                // need to mark as collision, not just close approach
                // keep going through data to accumulate objects in multiGp
                // and find appropriate tsFactor
                rollback = true;
            }

            // calc approach velocity - will need it
            // for either close approach or collision

            // criteria for timestep adjustment:
            // 1 - rate of momentum change
            //    calculate |dv|^2 / |v|^2 (unless v is close to 0), use as a criterion
            // 2 - rate of force angle change as a proportion of velocity
            //    calculate |dv|^2 - (dv dot v)^2 / |v|^2 = dv(normal)^2
            //    use dv(normal)^2 / |v|^2 as a criterion
            double tsFactorTmp = obj0->CheckTimestepCriteria();

            if (tsFactorTmp > tsFactor)
            {
                tsFactor = tsFactorTmp;
            }
        }

        // clean out the deltaVs for these objects (use collisions to remove paired data), 
        // set the timestep temporarily to mTStep / tsFactor
        // and run the calcs tsFactor times
        for (std::vector<CollisionData>::iterator cIt = collisions.begin();
            cIt != collisions.end();
            ++cIt)
        {
            Object* obj0 = cIt->mObj0;
            Object* obj1 = cIt->mObj1;
            obj0->RemoveDeltaVandInitialize(obj1);
            obj1->RemoveDeltaVandInitialize(obj0);
        }

        double tstep = mTStep / tsFactor;

        if (tsFactor > 1)
        {
            char msgbuf1[120];
            sprintf(msgbuf1, "\nProcessCollisions:  Temporarily set tstep to %g", tstep);
            LogSimpleMessage(msgbuf1);
        }

        // process all objects in this group using ministeps
        for (int i = 0; i < tsFactor; ++i)
        {
            CalcCollision(*mIt, tstep, tsFactor);
        }

        if (tsFactor > 1)
        {
            char msgbuf1[120];
            sprintf(msgbuf1, "\nProcessCollisions:   ... and tstep is back to %g", mTStep);
            LogSimpleMessage(msgbuf1);
        }
    }
}

void ObjectMgr::CalcStepForMultiObjects(std::vector<CollisionData>& multicolls, double tstep)
{
    // get set of objects and then process
    for (std::vector<CollisionData>::iterator mIt = multicolls.begin();
        mIt != multicolls.end();
        ++mIt)
    {
        Object* obj0 = mIt->mObj0;
        Object* obj1 = mIt->mObj1;
        double x0 = 0;
        double y0 = 0;
        double z0 = 0;
        obj0->GetLocation(x0, y0, z0);
        double x1 = 0;
        double y1 = 0;
        double z1 = 0;
        obj1->GetLocation(x1, y1, z1);
        double vx0 = 0;
        double vy0 = 0;
        double vz0 = 0;
        obj0->GetVelocity(vx0, vy0, vz0);
        double vx1 = 0;
        double vy1 = 0;
        double vz1 = 0;
        obj1->GetVelocity(vx1, vy1, vz1);
        double gmass0 = obj0->GetGMass();
        double gmass1 = obj1->GetGMass();

        // calculate delta vel based on half the acceleration,
        // add that to current vel, adjust loc based on vel, then add other half to vel
        // make sure to get delta vel in correct dir
        double x12 = x0 - x1;
        double y12 = y0 - y1;
        double z12 = z0 - z1;
        double rsq = (x12 * x12) + (y12 * y12) + (z12 * z12);
        double rdist = sqrt(rsq);

        // we want delta vels in 3 dirs
        // use dir12/r, so calc dv1 and dv2 with an extra r in the denom
        // calc 1/2 dv so we can get location change based on 1/2at**2
        double r3 = rdist * rsq;
        double dv1 = tstep * gmass1 / (2 * r3);
        double dv2 = tstep * gmass0 / (2 * r3);

        // register change in velocity with object
        // when pass is completed, have each object process the cumulative change
        obj0->AppendDeltaV(-x12 * dv1, -y12 * dv1, -z12 * dv1, obj1);
        obj1->AppendDeltaV(x12 * dv2, y12 * dv2, z12 * dv2, obj0);
    }
}

void FindCollisionGroups(std::vector<CollisionData>& collisions, std::vector<std::vector<CollisionData> >& multiCollisions)
{
    // objects may be in the initial list multiple times,
    // so need to find them and related ones and put them together in groups
    int numcollisions = collisions.size();
    std::set<Object*> processedObjs;

    // make sure we have all related close proximities grouped together
    for (int i = 0; i < numcollisions; ++i)
    {
        // if we've already processed an Object in this collision, 
        // both have been handled, so skip
        std::set<Object*>::iterator pIt = processedObjs.find(collisions[i].mObj0);

        if (pIt != processedObjs.end())
        {
            continue;
        }

        // set of Objects related to this collision data
        std::set<Object*> gpObjs;

        // grouped data
        std::vector<CollisionData> gp;
        gp.push_back(collisions[i]);

        // using a set so we can add indiscriminately and still have a unique collection
        gpObjs.insert(collisions[i].mObj0);
        gpObjs.insert(collisions[i].mObj1);

        // record so we don't process again
        processedObjs.insert(collisions[i].mObj0);
        processedObjs.insert(collisions[i].mObj1);

        // need a loop exit criterion:  use the size of collisions
        // when it doesn't change with a pass through the collision data, we're done
        int oldNum = 0;
        int newNum = 1;

        // keep checking the entire list as long as there is a new match
        // in practice, shouldn't amount to much
        while (oldNum != newNum)
        {
            oldNum = gpObjs.size();

            for (int j = 0; j < numcollisions; ++j)
            {
                if (i == j)
                {
                    continue;
                }

                // use tmp group, then add to gpObjs after loop
                std::set<Object*> objsToAdd;

                for (std::set<Object*>::iterator gIt = gpObjs.begin();
                    gIt != gpObjs.end();
                    ++gIt)
                {
                    if (collisions[j].ContainsObject(*gIt))
                    {
                        // add to group
                        objsToAdd.insert(collisions[j].mObj0);
                        objsToAdd.insert(collisions[j].mObj1);

                        // record so we don't process again
                        processedObjs.insert(collisions[j].mObj0);
                        processedObjs.insert(collisions[j].mObj1);

                        // add the collision to the group
                        gp.push_back(collisions[j]);
                    }
                }

                if (!objsToAdd.empty())
                {
                    for (std::set<Object*>::iterator aIt = objsToAdd.begin();
                        aIt != objsToAdd.end();
                        ++aIt)
                    {
                        gpObjs.insert(*aIt);
                    }
                }
            }

            newNum = gpObjs.size();
        }

        multiCollisions.push_back(gp);
    }
}

void ObjectMgr::CalcCollision(std::vector<CollisionData>& multicolls, double tstep, double tsfactor)
{
    // we're expected to have completed tstep and set v and loc for any colliding objects
    // but do nothing for others, unless we add them to list of collidedObjs
    double timeleft = tstep;
    std::set<Object*> objs;

    while (timeleft > 0)
    {
        double minCollisionTime = timeleft;
        bool hasCollision = false;

        for (std::vector<CollisionData>::iterator mIt = multicolls.begin();
            mIt != multicolls.end();
            ++mIt)
        {
            Object* obj0 = mIt->mObj0;
            Object* obj1 = mIt->mObj1;

            // record unique set of Objects to process at end
            objs.insert(obj0);
            objs.insert(obj1);

            // make sure mCurrent is initialized
            obj0->RemoveDeltaVandInitialize(obj1);
            obj1->RemoveDeltaVandInitialize(obj0);

            double proximityTime = 0;
            double dist = 0;
            bool isCollision = obj0->IntersectsInTimestep(obj1, proximityTime, dist, 1, timeleft);

            if (isCollision)
            {
                hasCollision = true;
                // calculate impact data:
                //   want ctr-ctr angle wrt each trajectory
                mIt->mDist = dist;
                mIt->mTimeToCollision = proximityTime;

                if (proximityTime < minCollisionTime)
                {
                    minCollisionTime = proximityTime;
                }
            }
        }

        // if we have some collisions, use the smallest time-to-collision as the next step
        // and then calc distances again
        if (hasCollision)
        {
            // take a step of length minCollisionTime
            // which puts at least 2 Objects in contact
            CalcStepForMultiObjects(multicolls, minCollisionTime);

            for (std::set<Object*>::iterator gIt = objs.begin();
                gIt != objs.end();
                ++gIt)
            {
                // and apply the calcs to location and velocity at each ministep
                // need the factor by which the original timestep is reduced, for applying deltavs
                assert((tstep - minCollisionTime) > 0);
                double factor = tsfactor * tstep / (tstep - minCollisionTime);
                (*gIt)->ProcessFractionalDeltaVs(objs, minCollisionTime, factor);
                (*gIt)->CommitTimeStepData();
            }

            timeleft -= minCollisionTime;

            // check to see if each pair is in contact
            // if so, fix velocities
            for (std::vector<CollisionData>::iterator mIt = multicolls.begin();
                mIt != multicolls.end();
                ++mIt)
            {
                if (mIt->mTimeToCollision == minCollisionTime)
                {
                    // collision:  reset velocities
                    Object* obj0 = mIt->mObj0;
                    Object* obj1 = mIt->mObj1;

                    // get ctr-ctr angle wrt each trajectory
                    // v0' = v0 - (2m1 / (m0 + m1)) [(v0 - v1) dot (x0 - x1)] (x0 - x1) / |(x0 - x1)|^2
                    // v1' = v1 + (2m0 / (m0 + m1)) [(v0 - v1) dot (x0 - x1)] (x0 - x1) / |(x0 - x1)|^2
                    // |(x1 - x0)| = r1 + r0
                    double x0 = 0;
                    double y0 = 0;
                    double z0 = 0;
                    obj0->GetLocation(x0, y0, z0);
                    double vx0 = 0;
                    double vy0 = 0;
                    double vz0 = 0;
                    obj0->GetVelocity(vx0, vy0, vz0);

                    double x1 = 0;
                    double y1 = 0;
                    double z1 = 0;
                    obj1->GetLocation(x1, y1, z1);
                    double vx1 = 0;
                    double vy1 = 0;
                    double vz1 = 0;
                    obj1->GetVelocity(vx1, vy1, vz1);

                    double dx = x0 - x1;
                    double dy = y0 - y1;
                    double dz = z0 - z1;

                    double dvx = vx0 - vx1;
                    double dvy = vy0 - vy1;
                    double dvz = vz0 - vz1;
                    double totmass = obj0->GetMass() + obj1->GetMass();
                    double radboth = obj0->GetRadius() + obj1->GetRadius();
                    double rad2 = radboth * radboth;
                    double m0factor = 2 * obj1->GetMass() / totmass;
                    double m1factor = 2 * obj0->GetMass() / totmass;
                    double dot = (dvx * dx) + (dvy * dy) + (dvz * dz);

                    // set new values reflecting collision in committed, and then reinitialize current
                    TimeStepData& tsd0 = obj0->GetCommitted();
                    tsd0.mVx = vx0 - m0factor * dot * dx / rad2;
                    tsd0.mVy = vy0 - m0factor * dot * dy / rad2;
                    tsd0.mVz = vz0 - m0factor * dot * dz / rad2;

                    TimeStepData& tsd1 = obj1->GetCommitted();
                    tsd1.mVx = vx1 + m1factor * dot * dx / rad2;
                    tsd1.mVy = vy1 + m1factor * dot * dy / rad2;
                    tsd1.mVz = vz1 + m1factor * dot * dz / rad2;
                    char msgbuf[120];
                    sprintf(msgbuf, "\nCalcCollision:  Objects %d and %d have collided",
                        obj0->GetIndex(), obj1->GetIndex());
                    LogSimpleMessage(msgbuf);
                }
            }
        }
        else
        {
            CalcStepForMultiObjects(multicolls, timeleft);

            for (std::set<Object*>::iterator gIt = objs.begin();
                gIt != objs.end();
                ++gIt)
            {
                // and apply the calcs to location and velocity at each ministep
                (*gIt)->ProcessFractionalDeltaVs(objs, timeleft, tsfactor);
                (*gIt)->CommitTimeStepData();
            }

            timeleft = 0;
        }
    }
}

void ObjectMgr::LogSimpleMessage(const char* msg)
{
    if (mLogFile.is_open())
    {
        mLogFile << msg;
    }

    printf(msg);
}
