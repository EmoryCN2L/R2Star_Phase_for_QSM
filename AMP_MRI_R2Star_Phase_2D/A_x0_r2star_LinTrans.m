classdef A_x0_r2star_LinTrans < LinTrans
	% A_x0_r2star_LinTrans:  Linear transform class for function handles combining step 3 and 4 

    % why is the Frobenius norm computed differently?
    % construct the Frobenius norm separately
	properties
		M	% output dimension: M=sm_num*echo_num*coil_num
		N	% input dimension: N=mat_sz(1)*mat_sz(2)*echo_num
		mat_sz	% the size of the 2D echo images
        wav_coef_len    % the wavelet coefficients length
		echo_num	% the number of coils
		echo_time	% the echo_time
		X0	% the X0 matrix
        R2star  % the R2star matrix
		E_R2star	% the exponential of R2star
        E_R2star_sq
        A   % the inverse wavelet transform
        Ah  % the wavelet transform
		S	% function handle for forward multiply-square
		St	% function handle for transpose multiply-square
        X0_num; % the number of nonzero entries in X0
		%FrobNorm    % 1/(M*N) times the squared Frobenius norm 
	end

	methods

		% Constructor
		function obj = A_x0_r2star_LinTrans(M,N,mat_sz,wav_coef_len,echo_num,echo_time,A,Ah,S,St)
			obj = obj@LinTrans;

			% manditory inputs
			if ~(isa(N,'numeric')&isa(M,'numeric'))
				error('First and second inputs must be integers')   
			end
			obj.M = M;
			obj.N = N;
			obj.mat_sz = mat_sz;
            obj.wav_coef_len = wav_coef_len;
			obj.echo_num = echo_num;
			obj.echo_time = echo_time;
            obj.A = A;
            obj.Ah = Ah;

			% optional inputs 
			if nargin > 8
				if isa(S,'double')&&(S>0)
					% 11th input "S" contains FrobNorm
					obj.FrobNorm = S;
				elseif (nargin > 8)&&(isa(S,'function_handle')&isa(St,'function_handle'))
					% 11th and 11th inputs are both function handles, S and St
					obj.S = S;
					obj.St = St;
				else
					error('Problem with the 9th & 10th inputs.  We need that either the fifth input is a positive number for FrobNorm, or that the fifth and sixth inputs are both function handles for S and St.')   
				end
			else

			end
			
			if (M~=mat_sz(1)*mat_sz(2)*echo_num)
				error('Dimension mismatch: M')
			end
			
			if (N~=wav_coef_len)
				error('Dimension mismatch: N')
			end
		end
		
		function obj = set.R2star(obj, R2star)
			obj.R2star = R2star;
            E_R2star = zeros([size(R2star,1) size(R2star,2) obj.echo_num]);
            for (i=1:obj.echo_num)
                E_R2star(:,:,i) = exp(-obj.echo_time(i)*R2star);
            end
            obj.E_R2star = E_R2star;
            obj.E_R2star_sq = E_R2star.^2;
		end
		
		function estFrob(obj)
            % this is only the Froubenius norm of one wavelet operator
			% approximate the squared Frobenius norm
			P = 2;      % increase for a better approximation
			obj.FrobNorm = 0;
			for p=1:P,   % use "for" since A may not support matrices 
				x_tmp = randn([obj.wav_coef_len 1]);
				norm_x_tmp = norm(x_tmp(:),'fro');
				y_tmp = obj.mult(x_tmp);
				obj.FrobNorm = obj.FrobNorm + (norm(y_tmp(:), 'fro')/norm_x_tmp).^2;
			end
			obj.FrobNorm = sqrt(obj.FrobNorm * (obj.N/P));
		end

		% Size
		function [m,n] = size(obj,dim)
			if nargin>1 % a specific dimension was requested
				if dim==1
					m=obj.M;
				elseif dim==2
					m=obj.N;
				elseif dim>2
					m=1; 
				else
					error('invalid dimension')
				end
			elseif nargout<2  % all dims in one output vector
				m=[obj.M,obj.N];
			else % individual outputs for the dimensions
				m = obj.M;
				n = obj.N;
			end
		end

		% Matrix multiply
		function y = mult(obj,x)

            x_im = obj.A(x);

			y=zeros([obj.mat_sz obj.echo_num]);
			for (j=1:obj.echo_num)
				y(:,:,j) = obj.E_R2star(:,:,j) .* x_im; 
			end
		end

		% Hermitian-transposed-Matrix multiply 
		function x = multTr(obj,y)
			for (j=1:obj.echo_num)
				y(:,:,j) = y(:,:,j) .* obj.E_R2star(:,:,j);
			end

            y_im = zeros(obj.mat_sz);
            for (j=1:obj.echo_num)
                y_im = y_im + y(:,:,j);
            end

            x = obj.Ah(y_im);
		end


		% Squared-Matrix multiply 
		function y = multSq(obj,x)
			if isempty(obj.FrobNorm)
				y = obj.S(x);
			else
                %y = ones([obj.mat_sz obj.echo_num])*((obj.FrobNorm^2/(obj.M*obj.N))*sum(x,'all'));
                y = (obj.FrobNorm^2/obj.M*x);
			end
		end


		% Squared-Hermitian-Transposed Matrix multiply 
		function x = multSqTr(obj,y)
			if isempty(obj.FrobNorm)
				x = obj.St(y);
			else
                %x = ones([obj.wav_coef_len 1])*((obj.FrobNorm^2/(obj.M*obj.N))*sum(y,'all'));
                x = (obj.FrobNorm^2/obj.N*y);
			end
		end

	end
end
