function [phi psi] = convertToPhiPsi(x,eta,u, g,alpha, H,F,sigma)
    u(isnan(u)) = 0;

    H_x     = eta - x*alpha;

    % u(x) => u(sigma)
    sigma_x = interp1(H, sigma, H_x);
    u_sigma = interp1(sigma_x, u, sigma);

    u_sigma(isnan(u_sigma))=0;

    % phi(x) => phi(sigma)
    phi = 2*g*eta;
    phi = interp1(sigma_x, phi, sigma);

    psi=F.*u_sigma;
    psi=psi(:);
end