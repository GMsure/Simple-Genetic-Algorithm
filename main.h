//
// Created by GMsure on 2021/4/20 0020.
//

#ifndef CLION_MAIN_H
#define CLION_MAIN_H
#include <iostream>
#include <atomic>
#include <ctime>

#include "Population.h"
#include "NatureSO.h"

//Cgt 最大无更新代数
#define CGT 75

//无变化精度
#define DLT 0.000001

//最佳适应度
#define PFIT 1


class Monitor{
public:
    Monitor();
    bool monitor(double VerFit);
    void setMFit(double fit);
    void update();
    double getPFit() const;
    void fitChange();

private:
    uint8_t cnt;//当前无变化计数
    double tmpFit;//前一次最佳适应度
    double mFit;
};




using namespace std;





#endif //CLION_MAIN_H
