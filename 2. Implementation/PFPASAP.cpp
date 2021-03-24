//
//  main.cpp
//  PFPASAP Project
//  Real-Tiem Embedded Systems Course - Winter 2020
//
//  Created by Mohammad Nasher on 13 April 2020.
//  Copyright © 2020 Mohammad Nasher. All rights reserved.
//  Contact: m[dot]nasher[at]iasbs[dot]ac[dot]ir

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cstdint>
#include <algorithm>
#include <filesystem>
#include <time.h>
#include <cmath>
#include <random>
#include <chrono>
#include <stdio.h>
#include <stdlib.h>

#define __DEBUG_COMMAND__ true

auto constexpr __DEF_INT_TEMP_VALUE__           = 2;
auto constexpr __POWER__FILE__ADDRESS__         = "Data.txt";
auto constexpr __DEF_BATTERY_VOLUME__           = 90;
auto constexpr __DEF_HYPER_PERIOD_SIZE__        = 60;
auto constexpr __DEF_MINIMUM_VALUE_BATTERY__    = 0;
auto constexpr __DEF_NUMBER_OF_PRIMES__         = 4;
//How Many Times You Want to Generate Tasksets for a Specific Utilization
auto constexpr __DEF_SIMULATION_DURATION__      = 1.0;

typedef long long int ll;

void _LogError(std::string _inp_Err){
    if(_inp_Err == ""){
        std::cout<<"<UKNOWN ERROR HAPPEND> PLEASE CHECKOUT THE REASON(S)...\n";
        return;
    }
    std::cout<<"[E]: <"<<_inp_Err<<"> PLEASE CHECKOUT THE REASON(S)...\n";
}

void _LogDetail(std::string _InpDetail){
    std::cout << "[D]: " << _InpDetail << std::endl;
}

class Task{
    private:
        //Worst-case execution time of task
        int _execTime;
        //Period of task
        int _periodTime;
        //Power consumption of task
        int _powerCons;
    
    public:
        //Task Constructor with no argument
        //  it will generate tasks with default parameters
        Task(){
            this->_execTime = __DEF_INT_TEMP_VALUE__;
            this->_periodTime = __DEF_INT_TEMP_VALUE__;
            this->_powerCons = __DEF_INT_TEMP_VALUE__;
        }
        
        //Task Constructor with param constructor
        //  it will generate tasks with given parameters
        Task(int _inp_exec, int _inp_period, int _inp_power){
            if(_inp_exec <= 0 || _inp_period <= 0 || _inp_power < 0){
                _LogError("INVALID ARGUMENT IN TASK::TASK(INT,INT,INT)");
                return;
            }
            this->_execTime = _inp_exec;
            this->_periodTime = _inp_period;
            this->_powerCons = _inp_power;
        }
        //Set execution time
        bool set_execTime(float _inp_value){
            if(_inp_value <= 0){
                _LogError("INVALID ARGUMENT IN TASK::SET_EXECTIME(FLOAT)");
                return false;
            }
            this->_execTime = _inp_value;
            return true;
        }
        //Get execution time
        int get_execTime(){
            return this->_execTime;
        }
        //Set period time
        bool set_periodTime(int _inp_value){
            if(_inp_value <= 0){
                _LogError("INVALID ARGUMENT IN TASK::SET_PERIODTIME(INT)");
                return false;
            }
            this->_periodTime = _inp_value;
            return true;
        }
        //Get period time
        int get_periodTime(){
            return this->_periodTime;
        }
        //Set power consumtion
        bool set_powerCons(int _inp_value){
            if(_inp_value <= 0){
                _LogError("INVALID ARGUMENT IN TAS::SET_POWERCONS(INT)");
                return false;
            }
            this->_powerCons = _inp_value;
            return true;
        }
        //Get power consumption
        int get_powerCons(){
            return this->_powerCons;
        }
};

class TaskSet{
private:
    double _Utilization;
    int _TaskSetSize;
    //A tasset is a set of tasks which specifies by 3 parameters -> T: (e, p ,pow)
    std::vector<Task> _Tasks;
    
    auto _randMake(int inp_size) {
        std::vector<double> result;
        std::random_device device;
        std::mt19937 generator(device());
        std::uniform_real_distribution<double> distribution(1,15);
        for (int i = 0; i < inp_size; i++){
            double rand_n =  distribution(generator);
            result.push_back(rand_n);
        }
        random_shuffle(result.begin(), result.end());
        return result;
    }

    void _genTasks(){
        for(int i = 0 ; i < this->_TaskSetSize; i++){
            Task* _temPTask = new Task();
            this->_Tasks.push_back(*_temPTask);
        }
        std::vector<double> _Priods;
        std::vector<double> _Uset;
        bool isValidUtilization = false;
        while(!isValidUtilization){
            //debug
            //std::cout << "create again Utilization " << std::endl;
            _Uset = this->_uUniFast();
            _Priods = this->_hyperPeriodLimitaion(__DEF_HYPER_PERIOD_SIZE__, this->_TaskSetSize);
            isValidUtilization = true;
            for(int j = 0 ; j < _Uset.size(); j++){
                if(_Uset[j] * _Priods[j] < 1)
                {
                    isValidUtilization = false;
                    break;
                }
            }
        }
        auto _randomPowers = this->_randMake(this->_TaskSetSize);
        for(int i = 0 ; i < this->_TaskSetSize; i++){
            this->_Tasks[i].set_powerCons(_randomPowers[i]);
           // std::cout << "U[i]: "<< _Uset[i] << "   period[i]: "<< _Priods[i] << "      exec[i]:" << _Uset[i]*_Priods[i] << std::endl;
            this->_Tasks[i].set_execTime(_Uset[i] * _Priods[i]);
            this->_Tasks[i].set_periodTime(_Priods[i]);
        }
        //Now all the task's parameter are ready
    }
    
    std::vector<double> _hyperPeriodLimitaionRand(int _Inp_Scale){
           //We produce some random numbers for now
           auto temp = this->_randMake(this->_TaskSetSize);
           for(auto i = 0; i < temp.size();i++)
               temp[i] = temp[i] * _Inp_Scale;
           return temp;
       }
    
    std::vector<int> _GetPrimesOf(int _Inp){
        if(_Inp <= 1){
            _LogError("INVALID ARGUMENT IN TASKSET::_GETPRIMESOF(INT)");
            return std::vector<int>(1,0);
        }
        std::vector<int> _Res;
        for(int i = 2 ; i < _Inp ; i++){
            bool isPrime = true;
            for(int j = 2 ; j < i; j++){
                if(i % j == 0)
                {
                    isPrime = false;
                    break;
                }
            }
            if(isPrime)
                _Res.push_back(i);
        }
        return _Res;
    }
    
    
    std::vector<double> _hyperPeriodLimitaion(int _Inp_UpperBound, int _Inp_Size){
        std::vector<double> _res;
        std::vector<std::vector<int>> _mMatrixVals;
        auto _Primes = this->_GetPrimesOf(_Inp_UpperBound);
        std::vector<int> _PrimesM;
        if(_Primes.size() > __DEF_NUMBER_OF_PRIMES__){
            for(int i = 0 ; i < __DEF_NUMBER_OF_PRIMES__; i++){
                _PrimesM.push_back(_Primes[i]);
            }
        }else{_LogError("LARGE NUMBER OF PRIMES HAS BEEN SET IN SIMULATION. IT CAN EFFECT OF RESULT OF SIMULATION");}
        //at the first step all of the cells are 1
        for(auto i = 0 ; i < __DEF_NUMBER_OF_PRIMES__; i++)
            _mMatrixVals.push_back(std::vector<int>(1,1));
        
        bool _CheckAgain = true;
        while(_CheckAgain){
            //std::cout << "while" << std::endl;
            _CheckAgain = false;
            //For each row in M Matrix
            for(auto i = 0 ; i < _mMatrixVals.size(); i++){
                //Create new Power and check with all other greatest powers in matrix
                int tmp = _PrimesM[i];
                for(auto j = 0 ; j < _mMatrixVals.size(); j++)
                    tmp = tmp * _mMatrixVals[j][_mMatrixVals[j].size() - 1];
                //One more power add to row
                if(tmp < _Inp_UpperBound){
                    _mMatrixVals[i].push_back(_mMatrixVals[i][_mMatrixVals[i].size() - 1] * _PrimesM[i]);
                    _CheckAgain = true;
                }
            }
        }
        //Now we have all matrix M Cell Values but what about accuracy
        std::vector<std::vector<int>> _mMatrix;
        for(auto i = 0 ; i < __DEF_NUMBER_OF_PRIMES__; i++)
            _mMatrix.push_back(std::vector<int>(1,1));
          //initialize Random Seed
          srand(time(NULL));
         //For each available Prime
        for(auto k = 0 ; k < _mMatrixVals.size(); k++){
            //For each Available power for ith prime
            for(int j : _mMatrixVals[k]){
                int sc = (rand() % _mMatrixVals[k].size()) + 1;
                for(int x = 0 ; x < sc; x++)
                    _mMatrix[k].push_back(j);
            }
        }

        int _tmpPeriod = 1;
        for(int i = 0 ; i < _Inp_Size; i++){
            _tmpPeriod = 1;
            //For each line of Matrix M
            for(std::vector<int> x : _mMatrix){
                auto seed = std::chrono::system_clock::now().time_since_epoch().count();
                std::uniform_int_distribution<int> distributionInteger(1, (int)x.size());
                std::default_random_engine generator(seed);
                int _Rand = distributionInteger(generator);
                _tmpPeriod = _tmpPeriod * x[_Rand - 1];
            }
            (_tmpPeriod > _Inp_UpperBound) ? _res.push_back(_Inp_UpperBound * 100):
                                             _res.push_back(_tmpPeriod * 100);
        }
        
        return _res;
    }
    
    std::vector<double> _uUniFast(){
        std::vector<double> UUniresult;
        auto random_list = this->_randMake(this->_TaskSetSize);
        std::random_device device;
        std::mt19937 generator(device());
        std::uniform_real_distribution<double> distribution(0,1);
        double sumU = this->_Utilization;
        double nextSumU;
        double rand_num;

        for (int j = 0; j < this->_TaskSetSize; j++) {
            rand_num = distribution(generator);
            nextSumU = sumU * pow(rand_num, (1.0 / (this->_TaskSetSize - j)));
            UUniresult.push_back(sumU - nextSumU);
            sumU = nextSumU;
        }
        UUniresult.push_back(sumU);
        return UUniresult;
    }
    
    int _ComputeLCM(int _Inp_First, int _Inp_Second){
           int _Max = (_Inp_First > _Inp_Second) ? _Inp_First : _Inp_Second;
           
          while (true){
              if (_Max % _Inp_First == 0 && _Max % _Inp_Second == 0)
                  break;
              else
                  ++_Max;
          }
           return _Max;
       }
    
public:
    //Generate taskset using utilization argument and
    //  number of tasks in the taskset
    TaskSet(float _inp_utilization, int _inp_size){
        if(_inp_utilization <= 0 || _inp_size <= 0){
            _LogError("INVALID ARGUMENT IN TASKSET::TASKSET(FLOAT,INT)");
            return;
        }
        this->_Utilization = _inp_utilization;
        this->_TaskSetSize = _inp_size;
        //Generate tasks with given parameters
        this->_genTasks();
    }
    
    //set and get functions
    float get_utilization(){
        return this->_Utilization;
    }
    bool set_utilization(double _inp_utilization){
        if(_inp_utilization <= 0){
            _LogError("INVALID ARGUMENT IN TASKSET::SET_UTILIZATION(INT)'");
            return false;
        }
        this->_Utilization = _inp_utilization;
        return true;
    }
    int get_size(){
        return this->_TaskSetSize;
    }
    bool set_size(int _inp_size){
        if(_inp_size <= 0 ||
           (this->_Tasks.size() != 0 && _inp_size != this->_Tasks.size())){
            _LogError("INVALID ARGUMENT IN TASKSET::SET_SIZE(INT)");
            return false;
        }
        this->_TaskSetSize = _inp_size;
        return true;
    }
    
    Task get_task(int _inp_index){
        if(_inp_index < 0 || _inp_index > this->get_size() - 1){
            _LogError("INDEX OUT OF RANGE IN TASKSET::GET_TASK(INT)");
            return Task();
        }
        return this->_Tasks[_inp_index];
    }
    
    int getRealHyperPeriod(){
        int _Res = 1;
        for(int i = 0 ; i < _Tasks.size(); i++){
            _Res = this->_ComputeLCM(_Res, this->_Tasks[i].get_periodTime());
        }
        return _Res;
    }
};


class SimulatorEngin{
private:
    std::string _FileAddress;
    TaskSet* _TastSet;
    std::vector<float> _PowerRate;
    std::vector<std::pair<double, int>> _SchedubilityRes;
    float _Battery;
    
    const std::vector<std::string> _explode(const std::string& _inp_string, const char& _inpSVal)
    {
        std::string _buffer{""};
        std::vector<std::string> _result;
        
        for(auto n: _inp_string)
        {
            if(n != _inpSVal) _buffer += n; else
            if(n == _inpSVal && _buffer != "") { _result.push_back(_buffer); _buffer = ""; }
        }
        if(_buffer != "") _result.push_back(_buffer);
        
        return _result;
    }

    bool _initPowerFormFile(){
        std::fstream _PowerFile;
        _PowerFile.open(this->_FileAddress);
        if(!_PowerFile){
            _LogError("POWER .TXT FILE NOT FOUND");
            return false;
        }
        std::string a;
        while(!_PowerFile.eof()){
            _PowerFile >> a;
            std::vector<std::string> v{_explode(a, ',')};
            if(v.size() == 0){
                return false;
            }
            try {
                int _rPower = std::stoi(v[v.size()-1]);
                if(_rPower < 0){
                    this->_PowerRate.push_back(0);
                }else{
                    // W/m*2 -> W/Cm*2 :)
                    this->_PowerRate.push_back(_rPower / 100);
                }
                    
            } catch (std::invalid_argument const &e) {
                _LogError("INVALID VALUES IN POWER DATASET FOUND BY SIMULATORENGINE::_INITPOWERFROMFILE(VOID)");
                continue;
            } catch(std::out_of_range const &e){
                _LogError("INTEGER OVERFLOW IN SIMULATORENGINE::_INITPOWERFROMFILE(VOID)");
            }
            
        }
        _PowerFile.close();
        if(this->_PowerRate.size() == 0){
            _LogError("NO VALUE IN POWER DATASET SIMULATORENGINE::_INITPOWERFROMFILE(VOID)");
            return false;
        }
        return true;
    }

    
    bool _IsValidAddress(const std::__fs::filesystem::path& p, std::__fs::filesystem::file_status s = std::__fs::filesystem::file_status{})
    {
        std::cout << p;
        if(std::__fs::filesystem::status_known(s) ? std::__fs::filesystem::exists(s) : std::__fs::filesystem::exists(p))
            return true;
        else
            return false;
    }
    
    bool _DoPFPasapWithBattery(TaskSet _inp_taskset){
        //in each time we can check ramaining job of task T by checking index T-1 of _RemainingJob
        std::vector<int> _RemainigJob;
        float _HPower = __DEF_BATTERY_VOLUME__;
        for(auto i = 0 ; i < _inp_taskset.get_size(); i++)
            _RemainigJob.push_back(_inp_taskset.get_task(i).get_execTime());
        
        for(auto _time = 1 ; _time < this->_PowerRate.size() * __DEF_HYPER_PERIOD_SIZE__; _time++){
            //at the begining of each hyper period
            if(_time % __DEF_HYPER_PERIOD_SIZE__ == 0){
                //Update the Power rate of hyper period
                _HPower = this->_PowerRate[_time / __DEF_HYPER_PERIOD_SIZE__];
            }
            //Check if there is any missed deadline
            for(auto _counter = 0 ; _counter < _RemainigJob.size(); _counter++){
                // because deadline = period ->
                // if time % deadline = 0 -> deadline missed cause there is a job that we
                //could not do it
                if(_time % _inp_taskset.get_task(_counter).get_periodTime() == 0){
                    // if there is any remaining job for task -> missed deadline
                    if(_RemainigJob[_counter] != 0){
                        std::cout << "No cause: " << _RemainigJob[_counter] <<"     "+ std::to_string(_time) <<std::endl;
                        return false;
                    }
                    //New Jobs Released -> we can do them until the next period
                    _RemainigJob[_counter] = _inp_taskset.get_task(_counter).get_execTime();
                }
                
            }
            //PFPasap Algorithm Start
               //find the higher priority available task at the moment
               int _AvailableIndex = -1;
               for(auto x = 0 ; x < _RemainigJob.size() ; x++)
                   if(_RemainigJob[x] != 0)
                       _AvailableIndex = x;
               if(_AvailableIndex == -1)
                   continue;
               //if have enough energy to do this job -> let's do it now :)
               if(this->_Battery + _HPower - __DEF_MINIMUM_VALUE_BATTERY__ >= _inp_taskset.get_task(_AvailableIndex).get_powerCons() / _inp_taskset.get_task(_AvailableIndex).get_execTime())
               {
                   //Do one job of task #x :D
                   _RemainigJob[_AvailableIndex]--;
                   //Manage Battery Power\
                   //If we need more than harvesting power we need to use (need - harvesting) power
                                                                                        //of battery
                   if(_inp_taskset.get_task(_AvailableIndex).get_powerCons()/ _inp_taskset.get_task(_AvailableIndex).get_execTime() >= _HPower){
                       this->_Battery -= (_inp_taskset.get_task(_AvailableIndex).get_powerCons()/ _inp_taskset.get_task(_AvailableIndex).get_execTime()) -  _HPower;
                   }
                   else{
                       //the task is gaining
                       //charge the battery if it's not fully charged
                       
                   }
               }
                //if do not have enough energy to do this job -> save the energy to do it in next moment
               else{
                   //if gain power in this moment battery will be fully charged
                   if(this->_Battery + _HPower >= __DEF_BATTERY_VOLUME__)
                       this->_Battery = __DEF_BATTERY_VOLUME__;
                   //else charge the battery till hpower end
                   else
                       this->_Battery += _HPower;
               }
            //PFPasap Algorithm End
        }
        return true;
    }
    
    bool _DoPFPasap(TaskSet _inp_taskset){
        //in each time we can check ramaining job of task T by checking index T-1 of _RemainingJob
        int _HyperPeriod = _inp_taskset.getRealHyperPeriod();
        #ifdef __DEBUG_COMMAND__
        std::cout << "HyperPriod = " << _HyperPeriod << std::endl;
        #endif
        std::vector<int> _RemainigJob;
        float _HPower = 0;
        this-> _Battery = (_HyperPeriod * 10) / 3;
        //std::cout << "Battery = " << this->_Battery << std::endl;
        for(auto i = 0 ; i < _inp_taskset.get_size(); i++)
            _RemainigJob.push_back(_inp_taskset.get_task(i).get_execTime());
        
        for(auto _time = 1 ; _time < this->_PowerRate.size() * _HyperPeriod; _time++){
            //at the begining of each hyper period
            if(_time % _HyperPeriod == 0){
                //Update the Power rate of hyper period
                _HPower = this->_PowerRate[_time / _HyperPeriod];
                //_HPower = 16;
            }
            //Check if there is any missed deadline
            for(auto _counter = 0 ; _counter < _RemainigJob.size(); _counter++){
                // because deadline = period ->
                // if time % deadline = 0 -> deadline missed cause there is a job that we
                //could not do it
                if(_time % _inp_taskset.get_task(_counter).get_periodTime() == 0){
                    // if there is any remaining job for task -> missed deadline
                    if(_RemainigJob[_counter] != 0){
                        #ifdef __DEBUG_COMMAND__
                        _LogDetail("DEADLINE MISSED BECAUSE " + std::to_string(_RemainigJob[_counter]) + " JOB OF " + std::to_string(_counter) + + "iTH TASK IN TASKSET AT TIME: "+std::to_string(_time));
                        #endif
                        return false;
                    }
                    //New Jobs Released -> we can do them until the next period
                    _RemainigJob[_counter] = _inp_taskset.get_task(_counter).get_execTime();
                }
                
            }
            //PFPasap Algorithm Start
               //find the higher priority available task at the moment
               int _AvailableIndex = -1;
               for(auto x = 0 ; x < _RemainigJob.size() ; x++)
                   if(_RemainigJob[x] != 0)
                       _AvailableIndex = x;
               if(_AvailableIndex == -1)
                   continue;
               //if have enough energy to do this job -> let's do it now :)
               if(this->_Battery + _HPower - __DEF_MINIMUM_VALUE_BATTERY__ >= _inp_taskset.get_task(_AvailableIndex).get_powerCons())
               {
                   //Do one job of task #x :D
                   _RemainigJob[_AvailableIndex]--;
                   //Manage Battery Power
                   this->_Battery = this->_Battery + _inp_taskset.get_task(_AvailableIndex).get_powerCons();
                   
               }
                //if do not have enough energy to do this job -> save the energy to do it in next moment
               else{
                    this->_Battery += _HPower;
               }
            //PFPasap Algorithm End
        }
        #ifdef __DEBUG_COMMAND__
        _LogDetail("TASKSET SUCCESSFULLY SCHEDULED");
        #endif
        return true;
    }
    
    
    //Log a TaskSet to a file
    void _Write2File(TaskSet _inp_ts, int index){
        std::ofstream _IOFile ("ts_" + std::to_string(_inp_ts.get_utilization()) + "_"
                              + std::to_string(index) + ".txt", std::ios::out | std::ios::app | std::ios::binary);
        for(auto i = 0 ; i < _inp_ts.get_size(); i++){
            _IOFile << _inp_ts.get_task(i).get_execTime() << "      ";
            _IOFile << _inp_ts.get_task(i).get_powerCons()<< "      ";
            _IOFile << _inp_ts.get_task(i).get_periodTime() << "        ";
            _IOFile << std::endl;
        }
    }
    
    void _PlotResult(){
        //Create and Write into Result file
        std::ofstream _IOFile ("PFPASAP_RES.txt", std::ios::out | std::ios::app | std::ios::binary);
        for(auto i = 0 ; i < this->_SchedubilityRes.size(); i++)
            _IOFile << this->_SchedubilityRes[i].first << "\t"<<this->_SchedubilityRes[i].second << std::endl;
        //Avoid Race Condition
        _IOFile.close();
        //Run Pyhton Plotter :)
        if (system(NULL)){
            system("python plot.py");
        }
        else
            _LogError("COMMAND PROCESSOR DOES NOT EXIST TO PLOT RESULT!");
        
    }

public:
    //Default Constructor With Void Argument
    SimulatorEngin(){
        this->_FileAddress = __POWER__FILE__ADDRESS__;
        this->_Battery = 0;
    }

    SimulatorEngin(std::string _inp_add){
        if(this->_IsValidAddress(_inp_add)){
            this->_FileAddress = _inp_add;
            return;
        }
        this->_FileAddress = __POWER__FILE__ADDRESS__;
        this->_Battery = 0;
    }
    
    void printres(){
        _LogDetail("RESULT OF SIMULATION: ");
        for(auto i = 0; i < this->_SchedubilityRes.size(); i++){
            _LogDetail(std::to_string((this->_SchedubilityRes[i].second / __DEF_SIMULATION_DURATION__) * 100.0) + "% SCHEDULABLE TASKSET FOR UTILIZATION = " + std::to_string(this->_SchedubilityRes[i].first));
        }
    }
    
    void start_engine(std::vector<float> _inp_utilization, int _inp_tsSize){
        if(_inp_utilization.size() == 0 || _inp_tsSize <= 0){
            _LogError("INVALID ARGUMANT IN SIMULATORENGINE::STARTENGINE(STD::VECTOR<FLOAT>,INT)");
            return;
        }
        if(!this->_initPowerFormFile())
            return;
        std::pair<double, int> resVal;
        for(auto i = 0 ; i < _inp_utilization.size(); i++){
            resVal.first = _inp_utilization[i];
            resVal.second = 0;
            this->_SchedubilityRes.push_back(resVal);
        }
         for(int  i = 0 ; i < _inp_utilization.size(); i++){
             if(_inp_utilization[i] <= 0 || _inp_utilization[i] > 1){
                 _LogError("INVALID ARGUMANT IN SIMULATORENGINE::STARTENGINE(STD::VECTOR<FLOAT>,INT)");
                 return;
             }
             for(auto j = 0 ; j < __DEF_SIMULATION_DURATION__ ; j ++){
                #ifdef __DEBUG_COMMAND__
                 _LogDetail("--------------- START TO SIMULATE  " + std::to_string(j) + "-th TASK OF UTILIZATION = " + std::to_string(_inp_utilization[i]));
                #endif
                 //Generate taskset with U = _inp_utilization[i]
                 this->_TastSet = new TaskSet(_inp_utilization[i], _inp_tsSize);
                 //Write Task set to a text file
                 this->_Write2File(*this->_TastSet, j + 1);
                 //Do PFPasap Scheduling for this task set
                 bool res = this->_DoPFPasap(*this->_TastSet);
                 if(res)
                     this->_SchedubilityRes[i].second++;
             }
         }
        //Print res
        this->printres();
        this->_PlotResult();
    }
};



int main(int argc, const char * argv[]) {
    
    
    SimulatorEngin* S = new SimulatorEngin();
    const clock_t begin_time = clock();
    S->start_engine(std::vector<float>{0.2, 0.4, 0.6, 0.8, 1.0} , 10);
    std::cout <<"SIMULATION DURATION IS: " <<float( clock() - begin_time ) /  CLOCKS_PER_SEC << " (s)" << std::endl;;
    
    return EXIT_SUCCESS;
}
