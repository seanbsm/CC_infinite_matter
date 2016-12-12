#ifndef MP_H
#define MP_H

class MP
{
public:
    MP();
    double  makeStateSpace  ();
    double  assym           (int p, int q, int r, int s);
    double  h0              (int p);
    double  f               (int p);
    int     kunique2        (int p, int q);
};

#endif // MP_H
