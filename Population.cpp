//
// Created by GMsure on 2021/4/23 0023.
//

#include "Population.h"


/// 必须使用该有参构造
/// \param sub_ptr 下标序列指针
/// \param num_ptr 数字对象序列指针
Population::Population(individual<int> * sub_ptr, individual<TYPE> * num_ptr) {
    fit = 0;
    size = sub_ptr->size();
    int size_tmp(num_ptr->size());
    if(size_tmp != size){
        std::cerr<<"Incorrect construct provided !"<<std::endl;
        std::cerr<<__FUNCSIG__<<std::endl<<"The size of the provided object must be consistent ! ( "<<size<<" != "<<size_tmp<<" )"<<std::endl;
    }
    if(size != PSIZE){
        std::cerr<<"Incorrect construct provided !"<<std::endl;
        std::cerr<<__FUNCSIG__<<std::endl<<"The content size of the supplied object must be "<<PSIZE<<std::endl;
    }

    individual<int> p1 = *sub_ptr;
    individual<TYPE> p2 = *num_ptr;
    for (int i = 0; i < size; i++) {
        C.push_back(mp2Code(p1.at(i)));
        EOS.push_back(p2.at(i));
    }


}

/// 得到基因
/// \param num 数字对象
/// \return 染色体基因
code Population::mp2Code(int num) {
    code m;
    for(int i =CODE_LEN - 1; i>=0; i--)
    {
        m[i] = (num>>i) & 1;//与1做位操作
    }
    return m;
}

/// 映射回数字sub
/// \param C 染色体基因
/// \return sub数字
int Population::mp2Num(code C) {
    int cArr = 0;
    for(int i = 0; i<CODE_LEN; i++){
        cArr += (C[i] << i);
    }
    return cArr;
}

/// 得到基因位置对应的下标
/// \param position 位置
/// \return 一个下标
int Population::getCodeSub(int position) {
    return mp2Num(C.at(position));
}

/// 更新适应度
/// \param fun 传入适应度函数重新计算的返回值
void Population::setFit(double fun) {
    fit = fun;
    if (fit == DBL_MAX){
        std::cerr<<"+INF"<<std::endl;
        fit=1;
    }
}

/// 得到当前适应度
/// \return
double Population::getFit() const {
    return  fit;
}

//打印种群信息
void Population::printSelfInfo() {
    std::cout<<"种群为:"<<std::endl;
    for (auto i: EOS){
        std::cout<<i<<" ";
    }
    std::cout<<std::endl;
    std::cout<<"适应度为:"<<fit<<std::endl;
    std::cout<<"总和为:"<<getSum()<<std::endl;

}

///  得到元素和
/// \return 当前种群元素和
TYPE Population::getSum() {
    TYPE sum = 0;
    for (auto i: EOS){
        sum += i;
    }
    return sum;
}

/// 得到位于position的基因
/// \param position 位置
/// \return 基因(sub)
code* Population::getCode(int position) {
    code * tmp (& C.at(position)) ;
    return tmp;
}

/// 得到当前种群位置对应的数字
/// \param position 位置
/// \return 一个下标
int Population::getCodeNum(int position) {
    return EOS.at(position);
}

/// 得到基因组
/// \return 种群基因组指针
individual<code> * Population::getGenome() {
    return &C;
}

/// 重建该种群
/// \param ptr 下标序列指针
/// \param num_ptr 数字对象序列指针
void Population::setNewPopulation(individual<int> *ptr, individual<TYPE> *num_ptr) {
    individual<int> p1 = *ptr;
    individual<TYPE> p2 = *num_ptr;
    for (int i = 0; i < size; i++) {
        C.push_back(mp2Code(p1.at(i)));
        EOS.push_back(p2.at(i));
    }
}

/// 重建该种群
/// \param ptr 基因序列指针
/// \param num_ptr 数字对象序列指针
void Population::setNewPopulation(individual<code> *ptr, individual<TYPE> *num_ptr) {
    individual<code> p1 = *ptr;
    individual<TYPE> p2 = *num_ptr;
    C.swap(p1);
    EOS.swap(p2);
//    for (int i = 0; i < size; i++) {
//        C.push_back(p1.at(i));
//        EOS.push_back(p2.at(i));
//    }
}




