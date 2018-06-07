#pragma once
#include <functional>
#include <vector>

class Barrier
{
    private:
        double contraction_start;
        double funnel_length;
        double tightest_length;
        double normal_width;
        double contracted_width;
    public:
     Barrier(double _contraction_start, double _funnel_length,
             double _tightest_length, double _normal_width,
             double _contracted_width) {
     
        contraction_start = _contraction_start;
        funnel_length = _funnel_length;
        tightest_length = _tightest_length;
        normal_width = _normal_width;
        contracted_width = _contracted_width;
     }
     double operator() (double x)
     {
         double tightest_start  = contraction_start + funnel_length;
         double tightest_end    = contraction_start + funnel_length + tightest_length;
         double contraction_end = contraction_start + 2*funnel_length + tightest_length;
         double width;
         if(x>contraction_start && x<tightest_start){
             width = (contracted_width - normal_width)/(funnel_length)*(x - contraction_start) + normal_width;
         }
         else if(x>=tightest_start && x<= tightest_end){
             width = contracted_width;
         }
         else if(x>tightest_end && x<contraction_end){
             width = (normal_width - contracted_width)/funnel_length *(x - tightest_end) + contracted_width;
         }
         else{
             width = normal_width;
         }

        return width;
     }
};

std::tuple<std::vector<double>, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>> solve_particle_model_simple(int, Barrier width);
std::tuple<std::vector<double>, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>> solve_particle_model(std::vector<double>&, std::vector<double>&, Barrier width);

