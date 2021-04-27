//
// Created by GMsure on 2021/4/23 0023.
//

#ifndef CLION_NATURESO_H
#define CLION_NATURESO_H
#include <random>

#include "Population.h"

//Zahl：初始数据集——全部初始个体数
#define ZAHL 50
//Pn：自然种群数量
#define PN 5
//Pc：交叉概率，一般取0.4 ~ 0.99；
#define PC 0.6789024
//Pm：变异概率，一般取0.0001~0.1。
#define PM 0.002

//去掉下面一行注释可减少调试输出
//#define DEBUG

//取绝对值
template<typename T>
inline T Abs(T x){
    return x < 0? -x: x;
}

inline std::random_device rd;
inline std::default_random_engine E{rd()};

[[maybe_unused]] inline std::bernoulli_distribution Pc(PC);
[[maybe_unused]] inline std::bernoulli_distribution Pm(PM);

inline bool SE(double x){
    std::bernoulli_distribution se(x);
    return se(E);
}

#ifndef DOUBLE
    inline std::uniform_int_distribution<TYPE> Num(1, 511);
#endif
#ifdef DOUBLE
    inline std::uniform_real_distribution<TYPE> Num(1.0, 511.0);
#endif
inline std::uniform_int_distribution<int> SUB(0, 7);
inline std::uniform_int_distribution<int> SubE(0, 19);


class NatureSO {
public:
    NatureSO();
    individual<TYPE>* specifyNum(int begin, int end);
    individual<int>* specifySub(int begin, int end);
    double funFit(TYPE sum) const;
    void updateFits(double f, int position);
    void updateAvg();
    TYPE getNum(int sub_);
    void save(individual<code>* GC);
    void select();
    void crossover(Population **P);
    TYPE getT10();
    void replace(Population **P);
    int maxFit();
    int minFit();

    individual<int> sub;//下标集

private:
    TYPE number[ZAHL];//初始对象集
    [[maybe_unused]] individual<code> tmpElite;//本轮精英种群
    double VerFit;//最佳适应度
    double AvgFit;//平均适应度
    double Fit[PN];//适应度集
    TYPE sumT10;//总值的1/10
    int par[2][PN];//父母标志

    bool duplicatePare(int a,int b);
    bool duplicateNum(TYPE a);
    static bool duplicatePopulate(individual<code> *p,int a);
    code cross2New(code a,code b);
    static int  mp2Num(code C);
    static bool mutation(code * GC);
    static bool REmutation(code * GC);


};


#endif //CLION_NATURESO_H
