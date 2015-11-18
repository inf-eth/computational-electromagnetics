function [ output_args ] = fxp( w, s0, delta )
%FXP Summary of this function goes here
%   Detailed explanation goes here
output_args = (fx(w+delta,s0)-fx(w,s0))/delta;
end

