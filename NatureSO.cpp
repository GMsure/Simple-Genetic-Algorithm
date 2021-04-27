//
// Created by GMsure on 2021/4/23 0023.
//

#include "NatureSO.h"

// 算子无参构造
NatureSO::NatureSO() {
    sumT10 = 0;
    TYPE tmp = 0;
    for(int i = 0; i<ZAHL; i++){
#ifndef DOUBLE
        tmp = Num(E) % 256;
#endif
#ifdef DOUBLE
        tmp = fmod(Num(E) , 256.0);
#endif
        if (!duplicateNum(tmp)){
            while (!duplicateNum(tmp))
                tmp = fmod(Num(E) , 256.0);
        }
        number[i] = tmp;
        sumT10 += number[i];
        sub.push_back(i);
    }// ZAHL 个 个体
    sumT10 /= 10;
    memset(par,0,sizeof (par));
    for(double & j : Fit){
        j = 0;
    }

}

///  得到[begin,end]该闭区间内自然对象集的全部个体
/// \param begin 起始位置
/// \param end  结束位置
/// \return 一个新的对象集合 -指针-
individual<TYPE>* NatureSO::specifyNum(int begin, int end) {
    if(begin >= end){
        std::cerr<<"The "<<__FUNCSIG__<<" function is being used incorrectly"<<std::endl;
        std::cerr<<"At Line_"<<__LINE__<<std::endl;
    }
    VerFit = 0;
   auto *tmp = new individual<TYPE>;
    for(int i = begin-1; i<end; i++){
        tmp->push_back(number[i]);
    }
    return tmp;
}

///  得到[begin,end]该闭区间内自然对象集的全部下标
/// \param begin 起始位置
/// \param end  结束位置
/// \return 一个新的对象集合 -指针-
individual<int> * NatureSO::specifySub(int begin, int end) {
    if(begin >= end){
        std::cerr<<"The "<<__FUNCSIG__<<" function is being used incorrectly"<<std::endl;
        std::cerr<<"At Line_"<<__LINE__<<std::endl;
    }
    VerFit = 0;
    auto *tmp = new individual<int>;
    for(int i = begin-1; i<end; i++){
        tmp->push_back(sub.at(i));
    }
    return tmp;
}

/// 适应度函数
/// 计算适应度
/// \param sum 元素和
/// \return
double NatureSO::funFit(TYPE sum) const {
    if (Abs<double>((sumT10 - sum*1.0)) <= 0.000001)
        return 1;
    else
        return (1.0 / pow(std::floor(Abs<double>((sumT10 - sum*1.0))) ,1.0));
}

///根据下标得到数字对象
/// \param sub_ 下标
/// \return 值
TYPE NatureSO::getNum(int sub_) {
    return number[sub_];
}

///
/// \param C
/// \return
int NatureSO::mp2Num(code C) {
    int cArr = 0;
    for(int i = 0; i<CODE_LEN; i++){
        cArr += (C[i] << i);
    }
    return cArr;
}

///
/// \param f
/// \param position
void NatureSO::updateFits(double f, int position) {
    if (f == +INFINITY){
        Fit[position] = 1;
#ifdef DEBUG
        std::cerr<<"+INF"<<std::endl;
#endif
    } else
        Fit[position] = f;
}

///得到最优适应度下标
/// \return 下标
int NatureSO::maxFit() {
    double tmp = Fit[0];
    int i = 0, j = 0;
    for (;i<PN; i++){
        if(Fit[i] > tmp){
            tmp = Fit[i];
            j = i;
        }
    }
    return j;
}

///得到最差适应度下标
/// \return 下标
int NatureSO::minFit() {
    double tmp = Fit[0];
    int i = 0, j = 0;
    for (;i<PN; i++){
        if(Fit[i] < tmp){
            tmp = Fit[i];
            j = i;
        }
    }
    return j;
}

///保存最优
/// \param GC 基因组指针
void NatureSO::save(individual<code>* GC) {
    tmpElite = *GC;
}

//更新平均适应度
void NatureSO::updateAvg() {
    double tmp = 0;
    for (auto i:Fit){
        tmp += i;
    }
    AvgFit = tmp/(PN*1.0);
}


/// 检查父母组合是否已存在
/// \param a 需要查重的元素组合 其一
/// \param b 需要查重的元素组合 其二
/// \return 重复 -0- 不重复 -1-
bool NatureSO::duplicatePare(int a, int b) {
    for(int i = 0; i<PN;i++){
        if((par[0][i] == a && par[1][i] == b) || (par[0][i] == b && par[1][i] == a)){
            return false;
        }
    }
    return true;
}

bool NatureSO::duplicatePopulate(individual<code> *p, int a) {
    for (auto i : *p){
        if (a == mp2Num(i)){
            return false;
        }
    }
    return true;
}

// 随机自然选择
// 会为下一代选择父母
void NatureSO::select() {
    int a=0,b=0;
    double prob = AvgFit;
    for (int i = 0; i< PN; i++){
        while (true){
            a = SUB(E) % PN;                         //生成下标1
            if(Fit[a] < AvgFit){                         //如果适应度小于均值
                prob = (Fit[a]*1.0) / AvgFit;       //则
                if (SE(prob)){                               //以一定的概率
                    a = SUB(E) % PN;                 //重新生成一次
                }
            }
            b = SUB(E) % PN;
            if(Fit[b] < AvgFit){
                prob = (Fit[b]*1.0) / AvgFit;
                if (SE(prob)){
                    b = SUB(E) % PN;
                }
            }
            if (a == b){
                continue;
            } else{
                if (duplicatePare(a,b)){
                    par[0][i] = a;
                    par[1][i] = b;
                    break;
                } else{
                    continue;
                }
            }
        }
    }
}

///
/// \param a
/// \param b
/// \return
code NatureSO::cross2New(code a, code b) {
    int tps = 0;
    double chance = 0;
    code tmp (b);
    for (int i = 1; i < CODE_LEN; ++i) {
        if(Pc(E)){
            tps = i;
        }
    }
    if (tps == 1){
        chance = number[(int)floor(Abs<double>(mp2Num(a) - mp2Num(b))) % ZAHL] - sumT10*1.0/PSIZE;
        if (Abs<double>(chance) <= 16){//在中位数上下16个数字之内再给一次机会
            for (int i = 0; i < CODE_LEN; ++i) {
                if(! Pc(E)){
                    tps = i;
                }
            }
        }
    }
    for (int i = tps; i < CODE_LEN; ++i) {
        tmp[i] = a[i];
    }
    return tmp;
}

///
/// \param P
void NatureSO::crossover(Population **P) {

    //individual<int> tmp_s;
    for (int i = 0; i < PN; ++i) {
        individual<code> tmp_c ;
        individual<TYPE> tmp_n;
        code M_c;
        int tmpa = 0, tmpb = 0;
        for (int j = 0; j < PSIZE; ++j) {
            M_c = cross2New(P[par[0][i]]->getGenome()->at(j),P[par[1][i]]->getGenome()->at(j));
            if(mutation(&M_c)){
                if(!duplicatePopulate(&tmp_c,mp2Num(M_c))){
                    while (!duplicatePopulate(&tmp_c,mp2Num(M_c)))
                        REmutation(&M_c);
                }
#ifdef DEBUG
                std::cout<<"M"<<std::endl;
#endif
            }
            if(!duplicatePopulate(&tmp_c,mp2Num(M_c))){
                while (!duplicatePopulate(&tmp_c,mp2Num(M_c)))
                    REmutation(&M_c);
            }
            tmp_c.push_back(M_c);
            tmp_n.push_back(getNum(mp2Num(tmp_c.at(j))%ZAHL));
            //tmp_s.push_back(mp2Num(tmp_c.at(j)));
        }
        if (Pc(E)){
            tmpa = SubE(E)%PSIZE;
            tmpb = SubE(E)%PSIZE;
            std::swap(tmp_c.at(tmpa),tmp_c.at(tmpb));
            std::swap(tmp_n.at(tmpa),tmp_n.at(tmpb));
        }

        P[i]->setNewPopulation(&tmp_c,&tmp_n);
        P[i]->setFit(funFit(P[i]->getSum()));
        updateFits(P[i]->getFit(),i);
        //P[i]->setNewPopulation(&tmp_s,&tmp_n);
    }
}

///
/// \param GC
/// \return
bool NatureSO::mutation(code* GC) {
    bool flag = false;
    for (int j = 0; j < CODE_LEN; ++j) {
        if (Pm(E)){
            GC->at(j) = ! GC->at(j);
            flag = true;
        }
    }
    return flag;
}
///
/// \param GC
/// \return
bool NatureSO::REmutation(code* GC) {
    bool flag = false;
    for (int j = 0; j < CODE_LEN; ++j) {
        if (Pc(E)){
            GC->at(j) = ! GC->at(j);
            flag = true;
        }
    }
    return flag;
}

///
/// \param P
void NatureSO::replace(Population ** P) {
    int i = minFit();
    individual<TYPE> tmp_i;
    for (int j = 0; j < PSIZE; ++j) {
        tmp_i.push_back(getNum(mp2Num(tmpElite.at(j))%ZAHL));
    }
    P[i]->setNewPopulation(&tmpElite,&tmp_i);
    P[i]->setFit(funFit(P[i]->getSum()));
    updateFits(P[i]->getFit(),i);
}

///
/// \return
TYPE NatureSO::getT10() {
    return sumT10;
}

///
/// \param a
/// \return
bool NatureSO::duplicateNum(TYPE a) {
    for (auto i : number) {
        if (a == i){
            return false;
        }
    }
    return true;
}
