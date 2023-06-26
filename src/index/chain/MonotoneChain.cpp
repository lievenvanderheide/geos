/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2006 Refractions Research Inc.
 * Copyright (C) 2001-2002 Vivid Solutions Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************
 *
 * Last port: index/chain/MonotoneChain.java rev. 1.15 (JTS-1.10)
 *
 **********************************************************************/

#include <geos/index/chain/MonotoneChain.h>
#include <geos/index/chain/MonotoneChainSelectAction.h>
#include <geos/index/chain/MonotoneChainOverlapAction.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/LineSegment.h>
#include <geos/geom/Envelope.h>
#include <geos/util.h>
#include <geos/util/Assert.h>
#include <geos/algorithm/Orientation.h>

#include <iomanip>

using namespace geos::geom;

namespace geos {
namespace index { // geos.index
namespace chain { // geos.index.chain

MonotoneChain::MonotoneChain(const geom::CoordinateSequence& newPts,
                             std::size_t nstart, std::size_t nend, void* nContext)
    : pts(&newPts)
    , context(nContext)
    , start(nstart)
    , end(nend)
    , env()
{}

const Envelope&
MonotoneChain::getEnvelope() const
{
    return getEnvelope(0.0);
}

const Envelope&
MonotoneChain::getEnvelope(double expansionDistance) const
{
    // This doesn't seem right. There's a single cached env, which is used with any value for 'expansionDistance'.
    if (env.isNull()) {
        env.init(pts->getAt<CoordinateXY>(start), pts->getAt<CoordinateXY>(end));
        if (expansionDistance > 0.0) {
            env.expandBy(expansionDistance);
        }
    }
    return env;
}

std::unique_ptr<CoordinateSequence>
MonotoneChain::getCoordinates() const
{
    return pts->clone();
}

void
MonotoneChain::select(const Envelope& searchEnv, MonotoneChainSelectAction& mcs) const
{
    computeSelect(searchEnv, start, end, mcs);
}

void
MonotoneChain::computeSelect(const Envelope& searchEnv,
                             std::size_t start0, std::size_t end0,
                             MonotoneChainSelectAction& mcs) const
{
    const CoordinateXY& p0 = pts->getAt<CoordinateXY>(start0);
    const CoordinateXY& p1 = pts->getAt<CoordinateXY>(end0);

    // terminating condition for the recursion
    if(end0 - start0 == 1) {
        mcs.select(*this, start0);
        return;
    }
    // nothing to do if the envelopes don't overlap
    if(!searchEnv.intersects(p0, p1)) {
        return;
    }
    // the chains overlap,so split each in half and iterate (binary search)
    std::size_t mid = (start0 + end0) / 2;

    // Assert: mid != start or end (since we checked above for end-start <= 1)
    // check terminating conditions before recursing
    if(start0 < mid) {
        computeSelect(searchEnv, start0, mid, mcs);
    }

    if(mid < end0) {
        computeSelect(searchEnv, mid, end0, mcs);
    }
}

/* public */
void
MonotoneChain::computeOverlaps(const MonotoneChain* mc,
                               MonotoneChainOverlapAction* mco) const
{
    computeOverlaps(start, end, *mc, mc->start, mc->end, 0.0, *mco);
}

/* public */
void
MonotoneChain::computeOverlaps(const MonotoneChain* mc, double overlapTolerance,
                               MonotoneChainOverlapAction* mco) const
{
    computeOverlaps(start, end, *mc, mc->start, mc->end, overlapTolerance, *mco);
}

/*private*/
void
MonotoneChain::computeOverlaps(std::size_t start0, std::size_t end0,
                               const MonotoneChain& mc,
                               std::size_t start1, std::size_t end1,
                               double overlapTolerance,
                               MonotoneChainOverlapAction& mco) const
{
    // terminating condition for the recursion
    if(end0 - start0 == 1 && end1 - start1 == 1) {
        mco.overlap(*this, start0, mc, start1);
        return;
    }

    // nothing to do if the envelopes of these subchains don't overlap
    if(!overlaps(start0, end0, mc, start1, end1, overlapTolerance)) {
        return;
    }

    // the chains overlap,so split each in half and iterate (binary search)
    std::size_t mid0 = (start0 + end0) / 2;
    std::size_t mid1 = (start1 + end1) / 2;

    // Assert: mid != start or end (since we checked above for
    // end-start <= 1)
    // check terminating conditions before recursing
    if(start0 < mid0) {
        if(start1 < mid1) {
            computeOverlaps(start0, mid0, mc, start1, mid1, overlapTolerance, mco);
        }
        if(mid1 < end1) {
            computeOverlaps(start0, mid0, mc, mid1, end1, overlapTolerance, mco);
        }
    }

    if(mid0 < end0) {
        if(start1 < mid1) {
            computeOverlaps(mid0, end0, mc, start1, mid1, overlapTolerance, mco);
        }
        if(mid1 < end1) {
            computeOverlaps(mid0, end0, mc, mid1, end1, overlapTolerance, mco);
        }
    }
}

/*private*/
bool
MonotoneChain::overlaps(const CoordinateXY& p1, const CoordinateXY& p2,
                        const CoordinateXY& q1, const CoordinateXY& q2,
                        double overlapTolerance)
{
    double maxq = std::max(q1.x, q2.x);
    double minp = std::min(p1.x, p2.x);
    if (minp > (maxq + overlapTolerance))
        return false;

    double minq = std::min(q1.x, q2.x);
    double maxp = std::max(p1.x, p2.x);
    if (maxp < (minq - overlapTolerance))
        return false;

    maxq = std::max(q1.y, q2.y);
    minp = std::min(p1.y, p2.y);
    if (minp > (maxq + overlapTolerance))
        return false;

    minq = std::min(q1.y, q2.y);
    maxp = std::max(p1.y, p2.y);
    if (maxp < (minq - overlapTolerance))
        return false;

    return true;
}

class MonotoneChain::BisectableIterator
{
public:
    BisectableIterator(const MonotoneChain& p_mc, bool forward);

    bool advance();

    bool bisect();

    const geom::CoordinateXY& getV0() const { return v0; }
    const geom::CoordinateXY& getV1() const { return v1; }

private:
    const MonotoneChain& mc;

    size_t index0;
    size_t index1;

    geom::CoordinateXY v0;
    geom::CoordinateXY v1;

    size_t stackHead;
    std::array<size_t, 64> stack;
};

MonotoneChain::BisectableIterator::BisectableIterator(const MonotoneChain& p_mc, bool forward)
    : mc(p_mc)
    , index0(forward ? mc.start : mc.end)
    , index1(forward ? mc.end : mc.start)
    , v0(mc.pts->getAt<geom::CoordinateXY>(index0))
    , v1(mc.pts->getAt<geom::CoordinateXY>(index1))
    , stackHead(0)
{
}

bool MonotoneChain::BisectableIterator::advance()
{
    if(stackHead == 0) {
        return false;
    }

    stackHead--;

    index0 = index1;
    index1 = stack[stackHead];

    v0 = v1;
    v1 = mc.pts->getAt<geom::CoordinateXY>(index1);

    return true;
}

bool MonotoneChain::BisectableIterator::bisect()
{
    if(index0 + 1 == index1 || index0 - 1 == index1) {
        return false;
    }

    stack[stackHead] = index1;
    stackHead++;

    index1 = (index0 + index1) / 2;
    v1 = mc.pts->getAt<geom::CoordinateXY>(index1);

    return true;
}

bool MonotoneChain::intersects(const MonotoneChain& mc) const
{
    geom::CoordinateXY left0 = pts->getAt<geom::CoordinateXY>(start);
    geom::CoordinateXY right0 = pts->getAt<geom::CoordinateXY>(end);
    bool chain0Forward = left0.compareTo(right0) == -1;
    if(!chain0Forward) {
        std::swap(left0, right0);
    }

    geom::CoordinateXY left1 = mc.pts->getAt<geom::CoordinateXY>(mc.start);
    geom::CoordinateXY right1 = mc.pts->getAt<geom::CoordinateXY>(mc.end);
    bool chain1Forward = left1.compareTo(right1) == -1;
    if(!chain1Forward) {
        std::swap(left1, right1);
    }

    int cmpResult = left0.compareTo(right1);
    if(cmpResult != -1) {
        return cmpResult == 0;
    }

    cmpResult = left1.compareTo(right0);
    if(cmpResult != -1) {
        return cmpResult == 0;
    }

    BisectableIterator it0(*this, chain0Forward);
    BisectableIterator it1(mc, chain1Forward);

    bool chain0BelowChain1;

    {
        bool advanceChain0 = left0.compareTo(left1) == -1;
        BisectableIterator& advancingIt = advanceChain0 ? it0 : it1;
        geom::CoordinateXY point = advanceChain0 ? left1 : left0;

        while(true) {
            if(point.compareTo(advancingIt.getV1()) == -1) {
                if(point.y < advancingIt.getV0().y && point.y < advancingIt.getV1().y) {
                    chain0BelowChain1 = !advanceChain0;
                    break;
                } else if(point.y > advancingIt.getV0().y && point.y > advancingIt.getV1().y) {
                    chain0BelowChain1 = advanceChain0;
                    break;
                } else if(!advancingIt.bisect()) {
                    int ordering = algorithm::Orientation::index(advancingIt.getV0(), advancingIt.getV1(), point);
                    if(ordering == algorithm::Orientation::COLLINEAR) {
                        return true;
                    }

                    chain0BelowChain1 = advanceChain0
                        ? (ordering == algorithm::Orientation::LEFT)
                        : (ordering == algorithm::Orientation::RIGHT);
                    break;
                }
            } else {
                bool advanceResult = advancingIt.advance();
                util::Assert::isTrue(advanceResult);
            }
        }
    }

    while(true)
    {
        bool advanceChain0 = it0.getV1().compareTo(it1.getV1()) == -1;
        BisectableIterator& advancingIt = advanceChain0 ? it0 : it1;
        BisectableIterator& otherIt = advanceChain0 ? it1 : it0;
        bool advancingBelowOther = advanceChain0 ? chain0BelowChain1 : !chain0BelowChain1;

        geom::CoordinateXY advancingV0 = advancingIt.getV0();
        geom::CoordinateXY advancingV1 = advancingIt.getV1();
        geom::CoordinateXY otherV0 = otherIt.getV0();
        geom::CoordinateXY otherV1 = otherIt.getV1();

        if(advancingBelowOther) {
            if(advancingV1.y > std::max(otherV0.y, otherV1.y)) {
                // advancingV1 is on the other side of the otherIt's current envelope, so there must be an intersection.
                return true;
            } else if(std::max(advancingV0.y, advancingV1.y) >= std::min(otherV0.y, otherV1.y)) {
                // The envelope of the advancing iterator and the other iterator intersect, so if possible, we should bisect,
                // otherwise, we've reached the level of a segment, so we directly compare against that.
                if(otherIt.bisect()) {
                    continue;
                } else if(algorithm::Orientation::index(otherV0, otherV1, advancingV1) != algorithm::Orientation::RIGHT) {
                    return true;
                }
            }
        } else {
            if(advancingV1.y < std::min(otherV0.y, otherV1.y)) {
                return true;
            } else if(std::min(advancingV0.y, advancingV1.y) <= std::max(otherV0.y, otherV1.y)) {
                if(otherIt.bisect()) {
                    continue;
                } else if(algorithm::Orientation::index(otherV0, otherV1, advancingV1) != algorithm::Orientation::LEFT) {
                    return true;
                }
            }
        }

        if(!advancingIt.advance()) {
            return false;
        }
    }
}

void MonotoneChain::printChain() const
{
    std::cout << std::setprecision(19) << std::endl;
    std::cout << "{" << std::endl;
    for(size_t i = start; i <= end; i++) {
        geom::CoordinateXY coords = pts->getAt<geom::CoordinateXY>(i);
        std::cout << "\t{" << coords.x << ", " << coords.y << "}," << std::endl;
    }
    std::cout << "};" << std::endl;
}

} // namespace geos.index.chain
} // namespace geos.index
} // namespace geos
