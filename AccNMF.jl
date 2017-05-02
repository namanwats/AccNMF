module AccNMF
export randinit , forb, MUW, MUH ,AccMU

# Function AccMU return factorised matrix W and H 
# using Accelated Multiplicatie Update
# randinit initialises matrix W and H
# forb is required to check the forbenius Norm
# MUH and MUW is used for multiplicative update of H and W

#####################################################
#		Random Initialization of W and H 			#
#####################################################
function randinit( X , k )
	p, n = size( X )
	T = eltype( X )

	W = rand( T , p , k )


	H = rand( T , k, n )
	return W , H
end

#####################################################
#				 Forbenius Norm 					#
#####################################################
function forb( E , M , N )
	sum = 0 ;
	for i = 1 : M
		for j = 1 : N
			sum += E[i,j] * E[i,j]
		end
	end
	e = sqrt( ( sum ) )
	return e
end


#####################################################
#			 Multiplicative Update					#
#####################################################


# %%%%%%%%%%%%%%% This function returns updated W %%%%%%%%%%%%%%%%% #

function MUW( X , W , A , B , R , M , N )


	temp2 = A ;
	temp3 = W * B ;


	Wp = W ;

	for i = 1 : M
		for j = 1 : R
			if ( isequal( temp2[i,j] , NaN) || isequal(temp2[i,j],Inf))
				temp2[i,j] = 1 ;
			end
			if ( isequal( temp3[i,j] , NaN ) || isequal( temp3[i,j] , Inf) )
				temp3[i,j] = 1 ;
			end
			if ( temp3[i,j] != 0 )
				W[i,j] = abs( W[i,j] * temp2[i,j] / temp3[i,j] ) ;
			else
				W[i,j] = 0 ;
			end
			if ( isequal( W[i,j] , NaN) || isequal( W[i,j] , Inf ) )
				W[i,j] = 1 ;
			end
		end
	end


	return W ;
end


# ************** This function returns updated W ***************** #

function MUH( X , A , B , H , R , M , N )
	Hp = H ;
	temp = A ;
	temp1 = B * H ;
	for i = 1 : R
		for j = 1 : N
			if ( isequal( temp[i,j] , NaN ) || isequal( temp[i,j] , Inf) )
				temp[i , j] = 1 ;
			end
			if(isequal(temp1[i,j] , NaN) || isequal(temp1[i,j] , Inf))
				temp1[i , j] = 1 ;
			end
			if(temp1[i,j] != 0)
				H[i,j] = abs(H[i,j] * (temp[i,j] / temp1[i,j]));
			end
			if(isequal(H[i,j] , NaN) || isequal(H[i,j] , Inf))
				H[i,j] = 1;
			end
		end
	end

	return H;
end


#########################################################
#					Calculate NMF 						#
#########################################################
function AccMU(X , R , alpha , Iter)
	println("NMF");
	
	W,H = randinit( X , R )	;	# Random initialization of W,H matrix


	Dims = size( X ) ;
	M = Dims[1] ;
	N = Dims[2] ;
	K = M * N ;

	pw = 1 + ( (K + N*R ) / ( M*R + M ) ) ;
	ph = 1 + ( (K + M*R ) / ( N*R + N ) ) ;


	LW = Int64(floor( 1 + alpha * pw )) ;	
	LH = Int64(floor( 1 + alpha * ph )) ;
	NW = zeros( M , R , LW ) ;
	NH = zeros( R , N , LH ) ;
	for k = 0 : Iter
		A = X * transpose(H) ;
		B = H * transpose(H) ;
		NW[: , : , 1] = W ;
		l = 1 ;
		T = eltype(W) ;
		for l = 2 : LW
			NW[: , : , l] = rand(T , M , R) ;
			NW[: , : , l] = MUW( X , NW[ : , : , l] , A , B , R , M , N) ;
			
			if(forb( NW[ : , : , l] - NW[ : , : , l - 1 ] , M , R ) <= 0.0001 * ( forb( NW[ : , : , 2] - NW[ : , : , 1 ] , M , R ) ) )
				break ;
			end
		end
		W = NW[ : , : , l ] ;


		#***** Calculating for H ****#

		A = transpose( W ) * X ;
		B = transpose( W ) * W ;
		NH[: , : , 1] = H;
		for l = 2 : LH
			NH[: , : , l] = rand(T , R , N ) ;
			NH[: , : , l] = MUH( X , A , B , NH[: , : , l] , R , M , N) ;
			if( forb( NH[: , : , l] - NH[ : , : , l - 1 ] , R , N ) <= 0.001 * (forb( NH[: , : , 2] - NH[: , : , 1] , R , N )  ) )
				break ;
			end
		end
		H = NH[: , : , l ] ;

	end
	
	return W,H;
end


end  #End of mudule


