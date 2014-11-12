

	%%%%%%%NANS 
	DUMMY_YY_NAN=find(isnan(DUMMY_YY(:)));  %%% removing NANS
	DUMMY_YY(DUMMY_YY_NAN)=[];
	DUMMY_XX(DUMMY_YY_NAN)=[];


	%%%%INFS%%%%%
	DUMMY_YY_inf=find(real(DUMMY_YY)==Inf);  % %%% removing infs
	DUMMY_YY(DUMMY_YY_inf)=[];
	DUMMY_XX(DUMMY_YY_inf)=[];
	

	%%%%
	
	DUMMY_YY_inf=find(real(DUMMY_YY)==-Inf);  % %%% removing infs
	DUMMY_YY(DUMMY_YY_inf)=[];
	DUMMY_XX(DUMMY_YY_inf)=[];


	%%% ZEROS %%%

	DUMMY_YY_zeros=find(real(DUMMY_YY)==0);   % %%% removing zeros
	DUMMY_YY(DUMMY_YY_zeros)=[];
	DUMMY_XX(DUMMY_YY_zeros)=[];
	
	%%%%%%


	%%%%%%  FINAL REAL NUMBERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	DUMMY_YY_real=real(DUMMY_YY);
	DUMMY_XX_real=real(DUMMY_XX);
	
	