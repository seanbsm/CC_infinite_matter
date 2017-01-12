#ifndef PM_H
#define PM_H

class PM
{
public:
    PM();
    double  makeStateSpace  ();
    double  assym           (int p, int q, int r, int s);
    double  assym_single    (int p, int q);
    double  h0              (int p);
    double  f               (int p);
    int     kUnique1        (int p);
    int     kUnique2        (int p, int q);
    int     kUnique3        (int p, int q, int r);
};

#endif // PM_H
