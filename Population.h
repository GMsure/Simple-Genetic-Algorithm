//
// Created by GMsure on 2021/4/23 0023.
//

#ifndef CLION_POPULATION_H
#define CLION_POPULATION_H
#include <map>
#include <vector>
#include <array>
#include <algorithm>
#include <iostream>

#define TYPE double
#define DOUBLE

//#define TYPE int


// 基因编码长度
#define CODE_LEN 6
// 种群个体数
#define PSIZE 10

using code = std::array<bool, CODE_LEN>;

template <typename T>
using individual = std::vector<T> ;

class Population {
public:
    explicit Population(individual<int> * ptr, individual<TYPE> * num_ptr);
    void setFit(double fun);
    double getFit() const;
    void printSelfInfo();
    TYPE getSum();
    code* getCode(int position);
    int getCodeSub(int position);
    int getCodeNum(int position);
    individual<code>* getGenome();
    void setNewPopulation(individual<int> * ptr, individual<TYPE> * num_ptr);
    void setNewPopulation(individual<code> * ptr, individual<TYPE> * num_ptr);


private:
    individual<code> C;//基因组 Vector
    double fit;//适应度
    uint32_t size;//种群大小
    individual<TYPE> EOS;//实体对象集Entity object

    static code mp2Code(int num);
    static int  mp2Num(code C);

};

#endif //CLION_POPULATION_H
