#include <iostream>
#include <math.h>
#include <cmath>
#include <eigen3/Eigen/Dense>

using namespace std;

template <typename T> //credit: http://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes#12399290
vector<size_t> sort_indexes(const vector<T> &v) {

    // initialize original index locations
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
         [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    return idx;
}

/*bool myfunction (int i, int j) {
  return (i==j);
}*/
bool contractor(int i, int j){ return i==j; }

int* returnInt(int i){
    return &i;
}

bool vecDelta(Eigen::VectorXi v1, Eigen::VectorXi v2){
    int dim1 = v1.rows();
    int dim2 = v2.rows();
    if (dim1 != dim2){
        cout << "dimensional error" << endl;
        return 0;
    }
    else{
        bool returnInt = 1;
        for (int i =0; i<dim1; i++){
            returnInt *= ( v1(i)==v2(i) );
        }
        return returnInt;
    }
}

int main()
{

    /*Eigen::MatrixXi states;
    Eigen::Matrix<int, 3, 10> m1;
    Eigen::Matrix<int, 3, 4> m2;*/

    Eigen::Vector4i v1(1,2,3,4);
    Eigen::Vector4i v2(1,2,3,4);

    //Eigen::VectorXi v = Eigen::VectorXi(1,2,3);

    Eigen::Matrix<int, 5, 1> vec;
    vec << 1,2,3,4,5;
    cout << vec.adjoint()*vec << endl;

    /*Eigen::VectorXi vec1;
    vec1 << 1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,5,5,5,5,5,5,5,5,6,6,6,6;
    Eigen::VectorXi vec2;
    vec2 << -1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,1,1,1,2,2,2,2,3,3,3,3,3,4,4,4,5,5,5,5,5,5,5,5,6,6,6,6,7,7,7;*/

    cout << vecDelta(v1+v1,v2+v2) << endl;

    /*for (int i= 0; i<1; i++){
        cout << i << endl;
    }*/

    /*m1 << 1,2,3,3,4,5,2,3,1,4,
            4,5,6,2,4,6,1,3,6,3,
            7,8,9,4,5,2,6,3,2,4;

    std::vector<int> v = {1,2,3,4,5,6,7,8};

    //Eigen::Vector4f v(1,2,3,4);
    v.erase( unique( v.begin(), v.end() ), v.end() );

    //cout << v.size() << endl;
    for (int i; i<v.size(); i++){
        cout << v[i] << endl;
    }
    int a = 5;
    cout << returnInt(a) << endl;

    m2 << 1,1,1,1,
            1,1,1,1,
            1,1,1,1;

    cout << ((m2.transpose())*m2).trace() << endl;*/
}
