#ifndef __SIMPLEMATCHERS_H_
#define __SIMPLEMATCHERS_H_

#include "model/planescan.h"
#include "match/matcher.h"

class EuklidMatcher : public Matcher
{
public:
    bool use_correspondence_min;
    EuklidMatcher(PlaneScan* ps, bool use_cm = false);
    void match();
};

class NormalMatcher : public Matcher
{
public:
    double eps_sim;
    NormalMatcher(PlaneScan* ps, double eps);
    void match();
};

#endif //__SIMPLEMATCHERS_H_