classdef A_2_echo_nw_LinTrans < LinTrans
	% FxnhandleLinTrans:  Linear transform class for function handles. 
    % no wavelet transform

	properties
		M	% output dimension: M=sm_num*echo_num*coil_num
		N	% input dimension: N=mat_sz(1)*mat_sz(2)*echo_num
		maps	% sensitivity maps, a 3d array with the third dimension being the coil indices
		maps_conj	% conjungate of maps
		mat_sz	% the size of the 2D echo images
		sm_num	% single channel measurement per echo
		echo_num	% the number of coils
		coil_num	% the number of echoes
		s_vect	% sampling vector
		A	% function handle for forward multiply
		Ah	% function handle for hermition-transpose multiply
		S	% function handle for forward multiply-square
		St	% function handle for transpose multiply-square
		%FrobNorm    % 1/(M*N) times the squared Frobenius norm 
	end

	methods

		% Constructor
		function obj = A_2_echo_nw_LinTrans(M,N,maps,mat_sz,sm_num,echo_num,coil_num,s_vect,A,Ah,S,St)
			obj = obj@LinTrans;

			% manditory inputs
			if ~(isa(N,'numeric')&isa(M,'numeric'))
				error('First and second inputs must be integers')   
			end
			obj.M = M;
			obj.N = N;
			obj.maps = maps;
			obj.maps_conj = conj(maps);
			obj.mat_sz = mat_sz;
			obj.sm_num = sm_num;
			obj.echo_num = echo_num;
			obj.coil_num = coil_num;
			obj.s_vect = s_vect;
			if ~(isa(A,'function_handle')&isa(Ah,'function_handle'))
				error('Third and fourth inputs must be function handles')   
			end
			obj.Ah = Ah;
			obj.A = A;

			% optional inputs 
			if nargin > 10
				if isa(S,'double')&&(S>0)
					% 11th input "S" contains FrobNorm
					obj.FrobNorm = S;
				elseif (nargin > 10)&&(isa(S,'function_handle')&isa(St,'function_handle'))
					% 13th and 14th inputs are both function handles, S and St
					obj.S = S;
					obj.St = St;
				else
					error('Problem with the 10th & 11th inputs.  We need that either the fifth input is a positive number for FrobNorm, or that the fifth and sixth inputs are both function handles for S and St.')   
				end
			else
				% approximate the squared Frobenius norm
				P = 10;      % increase for a better approximation
				obj.FrobNorm = 0;
				for p=1:P,   % use "for" since A may not support matrices 
				
					x_tmp = randn([mat_sz(1) mat_sz(2) echo_num]) + 1i * randn([mat_sz(1) mat_sz(2) echo_num]);
					norm_x_tmp = norm(x_tmp(:),'fro');
					y_tmp = obj.mult(x_tmp);
					obj.FrobNorm = obj.FrobNorm + (norm(y_tmp(:), 'fro')/norm_x_tmp).^2;
				end
				obj.FrobNorm = sqrt(obj.FrobNorm * (N/P));
			end
			
			if (M~=sm_num*echo_num*coil_num)
				error('Dimension mismatch: M')
			end
			
			if (N~=mat_sz(1)*mat_sz(2)*echo_num)
				error('Dimension mismatch: N')
			end
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
			y=zeros([obj.sm_num obj.echo_num obj.coil_num]);
			for (j=1:obj.echo_num)
                x_tmp = x(:,:,j);
				for (k=1:obj.coil_num)
					y(:,j,k) = obj.A(obj.maps(:,:,k).*x_tmp, obj.s_vect(:,j));
				end
			end
		end

		% Hermitian-transposed-Matrix multiply 
		function x = multTr(obj,y)
			x_tmp = zeros([obj.mat_sz obj.echo_num]);
			for (j=1:obj.echo_num)
				for (k=1:obj.coil_num)
					x_tmp(:,:,j) = x_tmp(:,:,j) + obj.mat_sz(1)*obj.mat_sz(2)*obj.maps_conj(:,:,k).*obj.Ah(y(:,j,k), obj.s_vect(:,j));
				end
			end
            x = x_tmp;
		end

		% Squared-Matrix multiply 
		function y = multSq(obj,x)
			if isempty(obj.FrobNorm)
				y = obj.S(x);
			else
				y = obj.echo_num*ones([obj.sm_num obj.echo_num obj.coil_num])*((obj.FrobNorm^2/(obj.M*obj.N))*sum(x,'all'));
			end
		end


		% Squared-Hermitian-Transposed Matrix multiply 
		function x = multSqTr(obj,y)
			if isempty(obj.FrobNorm)
				x = obj.St(y);
			else
				x = obj.echo_num*ones([obj.mat_sz(1) obj.mat_sz(2) obj.echo_num])*((obj.FrobNorm^2/(obj.M*obj.N))*sum(y,'all'));
			end
		end

	end
end
