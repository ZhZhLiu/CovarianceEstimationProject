function [ pdf ] = invGauPdf(x,lam,mu)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here

pdf = (lam/2/pi/x^3)^(1/2) * exp ( -lam * (x-mu)^2 / 2 / mu^2 / x );

end

