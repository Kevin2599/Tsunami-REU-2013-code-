%[sigma, F, H, H0, intF, dF, W, dW] = powerF(m,dsigma,sigmamax,ymax,g=1);
%
%Given a symmetric bey cross section defined by the power law f(y)=|y|^m,
%outputs the variety  of bayometrically determined variables required by
%our model.
%
%Inputs:
%-m is the power of the power function.
%-dsigma is the desired interval of the output sigma vector. Should not be
%smaller than 0.001.
%-sigmamax is the desired maximum of the output sigma vector.
%-ymax is a parameter that must be sufficiently large to allow for the
%desired sigmamax.
%-g is the value assigned to gravity. The default is 1.
%
%Outputs:
%-sigma is the independent output variable, with step size dsigma and
%maximum sigmamax.
%-F(sigma) is the smoothed desired function. It is calculated internally
%with a sigma step size of 0.001, smoothed, then, along with the remaining
%outputs, is cut down to use the sigma step size specified by dsigma. F is
%the stitched composite of two functions: close to the shore, it uses an
%analytic asymptotic expression for F, while farther out to sea it uses
%numeric methods. The stitching point is dynamically determined to be the
%point at which those two functions intersect.
%-H(sigma) is the total depth of water in the center of the channel, as a
%function of sigma.
%-H0(sigma) is the effective depth of water in the channel. It is
%calculated by dividing the area of the water in a cross-section by the
%width of the channel at the surface.
%-intF(sigma) is the integration of F(sigma). This is calculated using the
%identity intF = 2*g*H(sigma).
%-dF(sigma) is the smoothed first derivative of F(sigma). It is calculated using a
%second order central difference scheme, then smoothed using the matlab
%smooth() function. There is a jump at the stitching point between the
%asymptotic and numeric functions, which is smoothed using the fixit()
%function, which is a series of repeated linear smoothings. Thus, running
%trapF requires that the fixit.m function is accessible.
%-W(sigma) is a function of F and dF.
%-dW(sigma) is the derivative of W with respect to sigma.
%
%If you have questions, you can contact me (Jeremiah Harrington) at
%jeharrington@alaska.edu.

function [sigma,F,H,H0,intF,dF,W,dW] = powerF(m,dsigma,sigmamax,g)

%All of this can be solved explicitly.
%We start with an independent sigma.
sigma = [0:dsigma:sigmamax];
%F(sigma) is known for the power case.
F = (m/(m+1))*sigma;
%intF can be found by explicitly integrating this.
intF = (1/2)*(m/(m+1))*sigma.^2;
%So can dF.
dF = (m/(m+1))*ones(1,length(sigma));

%Note: intF=2gH.
H=intF/(2*g);
%F = 2sqrt(gH0).
H0 = (1/g)*(F/2).^2;

%Now, we find W and dW.
W = (2-dF)./F;
%Note, dF is constant.  Thus, dW = (2-dF)*(-1)*F^(-2)*dF
dW = (2-dF)*(-1).*(F.^(-2)).*dF;

end