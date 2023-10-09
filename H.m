function [val, error] = H(mu, w, E)%##########################################################################%Program written by M. Simcock, University of South Carolina, Columbia, SC%29208.  I can be reached at:  mn_sim@hotmail.com%%  The purpose of this program is to compute the value of Chandrasekhar's H%  function (sometimes called Ambartsumian's or Busbridge's H function)%  used in many diffuse reflectance and atmospheric/stellar light%  scattering and emission calculations.  The method for performing the%  calculation are used here comes from DWN Stibbs, RE Weir, Royal %  Royal Astronomical soc, vol 119, No5, 1959 page 512. The Function Invokes %  a  calculation using equation 21. One part of the method uses an %  infinite series; this program sums the first 100 terms of the series. %%  Inputs%  mu:      a value ranging from 0-1.  It is often the cosine of an angle.%  w:       the albedo, sometimes symbolized by the greek omega or the greek%           chi, the albedo also ranges from 0-1.  It is the ratio of the %           scattering cross-section to the sum of the scattering and %           absorption cross-sections.  %  E:       Specifies if an error estimate is to be output E = 0 (no error%           output) E~=0 and error estimate is output.%%  Outputs%  val:     This is the result of the H function%  error:   This is a estimate of the rounding error %%  Examples: %  [val, error] = H([0.7], [0.5]);%  [val, error] = H([0,0.05,0.10,1], [0.95,1]);%%  Notes:%  Regarding the "error" we estimate the rounding errors in the %  calculations by changing the value of the least significant figure in%  the input numbers by 1 (ie the sixteenth significant figure). The result%  is calculated and the difference taken to give an idea of the %  calculation error. %%  For some values of high mu and high w (eg mu = 0.96 and w = 0.9975) our %  calculated value is 2.603292264589538e+000 whereas the reference paper %  value is 2.603291.%  %  We suspect that the paper calculations were probably done using eight %  significant figures (the work was done in 1959) with the numbers reported %  to seven significant figures... so the slight disagreements in the paper %  values and the numbers calculated here may be due to their rounding errors %  rather than ours.% %  At the extremes where the function is valid (eg near w = 1) the calculation %  "error' jumps from 10e-16 (for moderate values) to 10e-8 (for more extreme %  values). This is estimated by changing the least significant figure in %  the number by 1. Users of this function should note that calculation errors %  become greater at the extremes of the function, from 1E-16 to 1E-8  (i.e%  when inputs close to zero and one are used).% %  A final point point regarding the "errors estimate" is that for some%  numbers, the rounding in the computation "fortuitously" compensates for the%  added error (i.e. it rounds to compensate the added error). As a result some%  errors come out as zero, even though there should be some rounding error%  somewhere, in any case teh rounding error is low in these cases and only %  affects teh least significant digits of the calculation output.% %##########################################################################format('long', 'e');%calculates H function with input values [val] = H_func(mu, w);error = 0;%calculates H function with values that have the least significant figure%altered by 1: this is done and an estimate of the propagation error that%occures in the calculation.if E ~= 0    mu_ = mu.*0;    for jj = 1:1:length(mu)        mu_(jj) = change_last_sig_dig(mu(jj));    end        w_ = w.*0;    for kk = 1:1:length(w)         w_(kk) = change_last_sig_dig(w(kk));    end         [val_] = H_func(mu_, w_);        error = abs(val - val_);end%##########################################################################function [val] = H_func(mu, w)    %Equation 21, integration to include zero point and pi/2.     tol = 1E-8; % quadrature tolerance    val = zeros(length(mu),length(w));    n_ = 6; %specifies the number of terms in eqn 14 to sum.             % we actually found that larger numbers of terms tended to            % increase the number and size of errors between Method 1 and            % Method 0.    mu_cut_off = (2.^0.5)-1; % his is the mu cut-off value defined in equation 15    for jj = 1:1:length(mu)    for kk = 1:1:length(w)             mu_ = mu(jj);          w_ = w(kk);               if mu_ == 0            val(jj,kk) = 1;        else                F = @(theta)I_func(w_,theta,mu_);            I_1 = (1./pi()).*quad(F,0,pi()./2,tol);                 if mu_ <= mu_cut_off                sum_part = 0;                for mm = 0:1:n_                    sum_part = sum_part + ((mu_.^((2.*mm)+1))./(((2*mm)+1).^2)); %summation part in equation 14.1                end                    I_2 = 0.5.*w_.*((0.5.*log(mu_).*log((1-mu_)/(1+mu_)))+ sum_part); % result of equation 14.1            else                sum_part = 0;                for mm = 0:1:n_                        sum_part = sum_part + ((1./(((2.*mm)+1).^2)).* (((1-mu_)./(1+(mu_))).^((2.*mm)+1))); %summation part of equation 14.2                end                I_2 = 0.5.*w_.*(((pi().^2)./8)-sum_part); %result of equation 14.2            end            if w_ == 1                I_4 = ((-2.*mu_)./pi()) .* (((3.*(1-w_))./w_).^0.5).*(((1-w_))).^0.5;               else                I_4 = ((-2.*mu_)./pi()) .* (((3.*(1-w_))./w_).^0.5).*...                atan((pi()./2).*((w_./(3.*(1-w_))).^0.5)); % equation 18            end                val(jj,kk) = exp(I_1 + I_2 + I_4);           end    end    end%##########################################################################function [F] = F_part(w_,theta,mu_)F = [];for kk = 1:1:length(theta)    if ((theta(kk) == 0)) && (w_ == 1) % condition of equation 15 in paper        F = [F,2.*mu_]; %equation 15, to compensate for theta -> 0 and w = 1, then F_part -> 2*mu    elseif theta(kk) == 0        F = 0;       else        F = [F,((w_.*(((theta(kk).*((1./sin(theta(kk))).^2))-(cot(theta(kk))))./...            (1-(w_.*theta(kk).*cot(theta(kk)))))).*(atan(mu_.*tan(theta(kk)))))];    endend%##########################################################################function [G] = G_part(w_,theta, mu_)G = ((pi()./2).*w_.*atan(mu_.*tan(theta)));%##########################################################################function [H] = H_part(w_,theta, mu_)H = [];for kk = 1:1:length(theta)    if ((theta(kk) == 0)) && (w_ == 1)        H = [H,0];        else        H = [H,((-2.*mu_.*(1-w_))./(1-w_+((1./3).*w_.*(theta(kk).^2))))];    endend%##########################################################################function [I] = I_func(w_,theta,mu_)I = F_part(w_,theta,mu_)... % f(theta), eq 11.2    -G_part(w_,theta, mu_)... % g(theta) eq 12    -H_part(w_,theta, mu_); %h(theta) equation 16 %##########################################################################function [san] = change_last_sig_dig(san)nas = num2str(san, '%.15e');    if str2num(nas(17)) == 9        nas(17) = '8';    else        if san == 1            nas = '0.9999999999999999';        else            nas(17)= num2str((str2num(nas(17)) + 1));        end    endsan = str2num(nas);