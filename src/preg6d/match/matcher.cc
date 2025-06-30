#include "match/matcher.h"

Matcher::Matcher(PlaneScan* p)
{
    setPlaneScan(p);
}

void Matcher::setPlaneScan(PlaneScan* p)
{
    this->ps = p;

}

