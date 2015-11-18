function [ output_args ] = jnup( nu, Z )
%JNUP Summary of this function goes here
%   Detailed explanation goes here
output_args = jnu(nu-1,Z)-((nu+1)/Z)*jnu(nu,Z);
end

