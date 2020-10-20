function [u, flag, relres, i] = glsqr(B, C, b, tol, maxit, iM, iN)
% function [U, FLAG, RELRES, ITER] = glsqr(B, C, b, TOL, MAXIT, iM, iN)
%
%   Attempts to solve 
%
%       C inv(N) B U = C inv(N) b
%
%   with a preconditioner defined by M (think M ~ C inv(N) B)
%
%   B and C are function handles that ideally compute (B*o) and (C*o)
%   iM & iN are function handles that ideally compute (M\o) and (N\o)
%
%   The solver exits whenever either holds:
%
%       a) MAXIT iterations performed
%       b) r1 = ||C N B U -  C N b||_{inv(M)} <= TOL
%
%   Diagnostic information:
%
%       FLAG    equals 0 if b) is satisfied, and 1 otherwise
%       RELRES  is NaN
%       ITER    is the number of iterations performed
%
%   Based on the LSQR algorithm of Paige and Saunders, 1982
%
%   The iterands (in particular b and U) may be cell arrays
%
%
%   See Sec 7.3 in
%
%       R. Andreev
%       Space-time discretization of the heat equation
%       Numerical Algorithms, 2014
%       (see README.txt for precise reference)

%   R. Andreev, 2012.10.18

	if (isempty(iM)); iM = @(u)(u); end
	if (isempty(iN)); iN = @(v)(v); end

	% Default values of the diagnostic output
	flag = 1;
	relres = NaN;

	% Initialization
	d = 0;
	[iNv, v, beta] = Normalize(b, iN);
	[iMw, w, alph] = Normalize(C(iNv), iM);
	rhoh = norm([alph beta]);
	u = my_times(0, iMw);
	delt = alph;
	gamm = beta;

	r1 = abs(delt) * gamm;
	if (r1 <= tol); flag = 0; return; end

	for i = 1:maxit
		d = my_plus(iMw, my_times(-alph*beta/rhoh^2, d));
		[iNv, v, beta] = Normalize(my_plus(B(iMw), my_times(-alph, v)), iN);
		[iMw, w, alph] = Normalize(my_plus(C(iNv), my_times(-beta, w)), iM);
		rhoh = norm([delt beta]);
		u = my_plus(u, my_times(+delt*gamm/rhoh^2, d));
		delt = -delt * alph / rhoh;
		gamm =  gamm * beta / rhoh;
		
		r1 = abs(delt) * gamm;
		if (r1 <= tol); flag = 0; return; end
		
		% % Note: r1 equals r0 from
		% [~,~,r0] = Normalize(my_plus(C(iN(B(u))), my_times(-1, C(iN(b)))), iM);
		% disp([r0 r1]);
	end
end

function [izSs, izs, z] = Normalize(s, iS)
	iSs = iS(s);
	z = sqrt(my_dot(s, iSs));
	izSs = my_times((1/z), iSs);
	izs = my_times((1/z), s);
end

%%% 
%   Provide custom functions for 
%      o) dot product
%      o) scalar multiplication
%      o) addition

function z = my_dot(X, Y)
	if iscell(X)
		assert(length(X) == length(Y), 'The cell arrays X and Y should be of the same length for cellwise dot product');
		z = 0;
		for i = 1:length(X)
			z = z + sum(sum(X{i} .* Y{i}));
		end
	else
		z = sum(sum(X .* Y));
	end
end

function Y = my_times(a, X)
	if iscell(X)
		for i = 1:length(X)
			Y{i} = a * X{i};
		end
	else
		Y = a * X;
	end
end

function Z = my_plus(X, Y)
	if iscell(X)
		assert(length(X) == length(Y), 'The cell arrays X and Y should be of the same length for cellwise addition');
		for i = 1:length(X)
			Z{i} = X{i} + Y{i};
		end
	else
		Z = X + Y;
	end
end
