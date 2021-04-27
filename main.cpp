//
// Created by GMsure on 2021/4/20 0020.
//

#include "main.h"

// 监视器无参构造
Monitor::Monitor() {
    cnt = 0;
    tmpFit = 0;
    mFit = 0;
}

/// 监视器
/// \param VerFit 上一次的适应度
/// \return 达到迭代目标 -1- 未达到精度或迭代要求 -0-
bool Monitor::monitor(double VerFit) {
    double tmp ;
    tmp = Abs<double>((VerFit - tmpFit));
    if (tmp < DLT){
        cnt ++;
    } else {
        cnt = 0;
    }
    if(cnt == CGT){
        return true;
    }
    return false;
}

//
void Monitor::fitChange() {
    if(tmpFit != mFit){
        cnt = 0;
    }
}

/// 设置当前最佳适应度
/// \param fit
void Monitor::setMFit(double fit) {
    mFit = fit;
}

//更新监视
void Monitor::update() {
    tmpFit =mFit;
    if (tmpFit == DBL_MAX){
        std::cerr<<"+INF"<<std::endl;
        tmpFit=1;
    }
}

/// 得到上一次的最佳适应度
/// \return
double Monitor::getPFit() const {
    return tmpFit;
}


int main(){

    system("chcp 65001");                                                        //指定终端编码UTF-8
    system("cls");
    system("cls");

    clock_t startTime,endTime;
    Monitor MNT;
    NatureSO NSO;
    Population *P[PN];
    long long int G = 0;

    for (int i = 1; i<ZAHL; i+=PSIZE){                                                               //初始化
        P[i/PSIZE] = new Population(NSO.specifySub(i+0, i+9), NSO.specifyNum(i+0, i+9));
        P[i/PSIZE]->setFit(NSO.funFit(P[i/PSIZE]->getSum()));
        NSO.updateFits(P[i/PSIZE]->getFit(),i/PSIZE);
        cout<<"第"<<i/PSIZE+1;
        P[i/PSIZE]->printSelfInfo();
    }

    cout<<"以上为初始种群，请等待算法求解"<<endl;
    cout<<"以上为初始种群，请等待算法求解"<<endl;
    cout<<"以上为初始种群，请等待算法求解"<<endl;
    startTime = clock();
    MNT.setMFit(P[NSO.maxFit()]->getFit());
    NSO.updateAvg();
    while (!MNT.monitor(MNT.getPFit())){

        NSO.save(P[NSO.maxFit()]->getGenome());                                     //精英保留
        NSO.select();                                                                                       //选择
        NSO.crossover(P);                                                                              //交叉及变异
        NSO.replace(P);                                                                                  //替换
        MNT.setMFit(P[NSO.maxFit()]->getFit());
        MNT.fitChange();
        MNT.update();                                                                                    //更新

#ifdef DEBUG
        for (int i = 0; i < PN; ++i) {
            cout<<"第"<<i+1;
            P[i]->printSelfInfo();
        }
#endif
        G++;
        NSO.updateAvg();
    }

    endTime = clock();
    cout<<"迭代已结束，用时"<<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s\n" <<"最优解集为:"<<endl;
    P[NSO.maxFit()]->printSelfInfo();
    cout<<"总和的十分之一为："<<NSO.getT10()<<endl;
    cout<<"共迭代了"<<G<<"代达到此效果"<<endl;


    return 0;
}
