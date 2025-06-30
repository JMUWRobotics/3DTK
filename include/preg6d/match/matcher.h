#ifndef __MATCHER_H_
#define __MATCHER_H_

class PlaneScan; // forward-declaration of planescan class

class Matcher
{

protected:
    PlaneScan* ps;

public:
    Matcher(PlaneScan* p);
    virtual ~Matcher() {};

    void setPlaneScan(PlaneScan* p);

    virtual void match() = 0;
};

#endif // __MATCHER_H_