% [Phi_nm1] = convertInitEta(eta,g)
%
% This program takes a physical eta function at t=0, and converts it into
% the abstract Phi at lambda=0 to be used as initial conditions in the
% model. This assumes that u=0 for all x at t=0.
%
% Note that the relation used in this program, that Phi_nm1 = 2*g*eta is
% derived from the back substitution
% eta = (phi - u^2)/(2*g)
% by assuming that at the initial time u = 0 and solving for eta.
%
% Inputs:
% eta - A vector based on x that defines the perturbation of the water
% surface at the initial time.
% g - The value assigned to gravity in our model.
%
% Outputs:
% Phi_nm1 - A vector across sigma that is an initial condition for the
% differential equation.

function [ Phi_nm1 ] = convertInit( eta,g )
Phi_nm1 = 2*g*eta;
end