#pragma once
#include<iostream>
#include<string>
#include<vector>


void subset(std::vector<int> v, int size, int left, int index, std::vector<int>& l,std::vector<std::vector<int>>& subsets){
    
    if(left==0){
        subsets.push_back(l);
        return;
    }
    for(int i=index; i<size;i++){
        l.push_back(v[i]);
        subset(v,size,left-1,i+1,l,subsets);
        l.pop_back();
    }

}     
