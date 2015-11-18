function [ output_args ] = h2nup( nu, Z )
%H2NUP Summary of this function goes here
%   Detailed explanation goes here
output_args = h2nu(nu-1,Z)-((nu+1)/Z)*h2nu(nu,Z);
end

