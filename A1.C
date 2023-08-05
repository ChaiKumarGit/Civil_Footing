
/*     DETEMINATION OF BEARING CAPACITY AND SETTLEMENTS      */

#include<stdio.h>
#include<math.h>
#include<conio.h>


main()
{

float fie,Nc,Nq,Nr,No,Sc,Sq,Sr,Dc,Dq,Dr,Ic,Iq,Ir,C,y,Df,b,L,A,k,Wl,H,V,Rq,Ry,Zq,Zy,Qult,Qs,FS,S1=0,S2,S12=100,Cc,e,Po,P,Wel,cen,TFM=0,Q,Qs1,FM[30],Cs,Cr,Is,Es,m,TIs,t1,t,y1,H1,S22,cen1,Val;
int a,op,i,nl,ch,jump=0,Bch,N,bstart=0,sstart=0;

FILE *fp,*fq;

fp=fopen("input1.txt","r");
fq=fopen("output1.txt","w");
fscanf(fp,"%f%f%f%d%f%f%d%d%f%f%f%f%f%f%d%f%f%f%d%f%f%d%f%f",&y,&Df,&Wl,&a,&b,&L,&Bch,&op,&C,&fie,&H,&V,&FS,&Q,&ch,&Cc,&e,&t,&nl,&m,&Es,&N,&t1,&y1);


       
 SF:
	



		if(a==1)
		{
	      A=b;
       	 Is=1.26;
		}
        else if(a==2)
		{
 		 A=(0.7854*b*b);
		 Is=0.73;
		}
        else if(a==3)
		{
          A=(b*L);
		 TIs=L/b;

		if(TIs==1)
			Is=0.82;
		else if(TIs==2)
			Is=1.0;
		else if(TIs==5)
			Is=1.22;
		else if(TIs==10)
			Is=1.26;
		else
			Is=1.0;

		}

        else if(a==4)
		{
         A=(b*b);
	   	 Is=0.82;
		}
        else
		{
 		goto SF;
		}

		
		if(Bch!=1)
		{
			jump=1;
			goto DSF;

		}
		else
		{

	     	if(Wl>(Df+b))
			{
			    Rq=1.0;
				Ry=1.0;
 			}
		    else if(Wl>=Df&&Wl<=(Df+b))
			{
			    Zy=(Wl-Df);
			    Rq=1;
	            Ry=0.5*(1+(Zy/b));
			}
	    	else
			{
			    Zq=Wl;
		        Rq=0.5*(1+(Zq/Df));
			    Ry=0.5;
			}

             if(fie==0)
				Nc=5.14 , Nq=1.00 , Nr=0.00  ;
			else
			{
                No=(tan((45+(fie/2))*(3.142/180)))*(tan((45+(fie/2))*(3.142/180)));
     		    Nq=pow(2.718,(3.142*tan(fie*(3.142/180))))*No;
			    Nc=(Nq-1)*(1/tan(fie*(3.142/180)));
			    Nr=1.8*(Nq-1)*tan(fie*(3.142/180));
			}

		}
		
	switch(op)
	{
	case 1:
 
        if(a==1)
	    	Qult=(C*Nc)+( y*Df*Nq*Rq)+(0.5*y*b*Nr*Ry);
        else if(a==2)
            Qult= (1.3*C*Nc)+( y*Df*Nq*Rq)+(0.3*y*b*Nr*Ry);
	    else if(a==3)
	    	Qult= (C*Nc)+( y*Df*Nq*Rq)+(0.5*y*b*Nr*Ry);
        else
            Qult= (1.3*C*Nc)+( y*Df*Nq*Rq)+(0.4*y*b*Nr*Ry);


          
 			Qs=(Qult/FS)+((y*Df*(FS-1))/FS);

    		

		break;

	case 2:


  if(a==1)
  {
	   Sc=1.0; Sq=1.0; Sr=1.0;
  }
    else if(a==2)
	{
	   Sc=1.3; Sq=1.2; Sr=0.6;
 	}
    else if(a==3)
	{
       Sc=(1+(0.2*(b/L)));
       Sq=(1+(0.2*(b/L)));
	   Sr=(1-(0.4*(b/L)));

	}
    else
	{
	   Sc=1.3; Sq=1.2; Sr=0.8;
	}


    Dc=(1+(0.35*(Df/b)));
	Dq=Dc;
	Dr=1.0;
    
 	Ic=(1-(H/(V+A*C*(1/tan((fie)*(3.142/180))))))*(1-(H/(V+A*C*(1/tan((fie)*(3.142/180))))));
	Iq=Ic;
	Ir=(Iq*Iq);

                 Qult= (C*Nc*Sc*Dc*Ic)+( y*Df*Nq*Sq*Dq*Iq*Rq)+(0.5*y*b*Nr*Sr*Dr*Ir*Ry);

	
     	Qs=(Qult/FS)+((y*Df*(FS-1))/FS);
     
  break;
	}
    fprintf(fq,"\n\n WIDTH OF THE FOOTING IS %f", b);
	fprintf(fq,"\n\tUltimate Bearing capacity(in kPa): %f",Qult);
			
	fprintf(fq,"\n\tSafe Bearing capacity(in kPa): %f",Qs);
			

DSF:
	if(jump==1)
	{
 	}

bbegin:
		if(sstart==0)
		{
	
		if(a==1)
		{
	      A=b;
       	 Is=1.26;
		}
        else if(a==2)
		{
 		 A=(0.7854*b*b);
		 Is=0.73;
		}
        else if(a==3)
		{
          A=(b*L);
		}
		else if(a==4)
		{
         A=(b*b);
	   	 Is=0.82;
		}
        else
		{
 		goto SF;
		}

		}
	  if(sstart!=0)
		{

		if(a==1)
		{
	      A=b;
       	 Is=1.26;
		}
        else if(a==2)
		{
 		 A=(0.7854*b*b);
		 Is=0.73;
		}
        else if(a==3)
		{
          A=(b*L);
		}
		else if(a==4)
		{
         A=(b*b);
	   	 Is=0.82;
		}
        else
		{
 		goto SF;
		}

		}



DL:

	switch(ch)
	{
		case 1:

	 
            
		    Qs1=Q/A;
	        Wel=t/nl;
	        cen=Wel/2;

			for(i=1;i<=nl;i++)
			{
				if((cen+Df)<=Wl)
				{
					Po=(cen+Df)*y;
		            P=Q/((b+cen)*(b+cen));
		            FM[i]=log10((Po+P)/Po);
		            TFM=TFM+FM[i];
		            cen=cen+Wel;
				}
	            else
				{
		            Po=((y-9.81)*(Df+cen-Wl))+(y*Wl);
	            	P=Q/((b+cen)*(b+cen));
					FM[i]=log10((Po+P)/Po);
	            	TFM=TFM+FM[i];
	            	cen=cen+Wel;
				}

			}
            

	        S1=1000*((Cc*Wel)/(1+e))*TFM;
			S2=1000*Qs1*b*((1-(m*m))/Es)*Is;
			S12=S1+S2;
			
			TFM=0;
			if(bstart==0)
			{
				fprintf(fq,"\n\nProposed Load on the Footing(in kN): %f",Q);
				fprintf(fq,"\n\n\tCONSOLIDATION SETTLEMENT OF CLAY(in mm)= %f",S1);
				fprintf(fq,"\n\n\tIMMEDIATE SETTLEMENT OF CLAY(in mm)= %f\n",S2);
		        fprintf(fq,"\n\n\tTOTAL SETTLEMENT OF CLAY(in mm)= %f\n",S12);
                fprintf(fq,"\n\tAllowable Bearing Pressure(in kPa): %f",Qs1);
			}
		
			if(S12>=75)
			{
				bstart=1;
				if(bstart==1)
				b=b+0.05;
				goto bbegin;
				
			}
			if(bstart==1)
			{
               	fprintf(fq,"\n\n\n**************************************************************");
				fprintf(fq,"\n\n\nAS THE SETTLEMENT HAS TO BE RESTRICTED TO 75mm AS PER IS-CODES");
		    	fprintf(fq,"\n\n\nWidth of the footing is changed to %f", b);
	            fprintf(fq,"\n\n\tCONSOLIDATION SETTLEMENT OF CLAY(in mm)= %f",S1);
				fprintf(fq,"\n\n\tIMMEDIATE SETTLEMENT OF CLAY(in mm)= %f\n",S2);        
		        fprintf(fq,"\n\n\tTOTAL SETTLEMENT OF CLAY(in mm)= %f\n",S12);
             	fprintf(fq,"\n\tAllowable Bearing Pressure(in kPa): %f",Qs1);
				}

			break;

		case 2:

			Qs1=Q/A;
		    H1=(b*3);
		    cen1=(H1/2);
			if((cen1+Df)<=Wl)
			{
				Po=(cen1+Df)*y;
		        P=Q/((b+cen1)*(b+cen1));
		        Val=log10((Po+P)/Po);
			}
	        else
			{
		        Po=((y-9.81)*(Df+cen1-Wl))+(y*Wl);
		        P=Q/((b+cen1)*(b+cen1));
		        Val=log10((Po+P)/Po);
			}

			Cr=400*N;
		    Cs=((1.5*Cr)/Po);
		    S22=1000*(H1/Cs)*Val;

		
			if(bstart==0)
			{
				fprintf(fq,"\n\nProposed Load on the Footing(in kN): %f",Q);
	            fprintf(fq,"\n\n\tIMMEDIATE SETTLEMENT OF SAND(in mm)= %f\n",S22);
			    fprintf(fq,"\n\tAllowable Bearing Pressure(in kPa): %f",Qs1);
			}
		
			if(S22>=75)
			{
				bstart=1;
				if(bstart==1)
				b=b+0.05;
				goto bbegin;
				
			}
			if(bstart==1)
			{
				fprintf(fq,"\n\n\n**************************************************************");
				fprintf(fq,"\n\n\nAS THE SETTLEMENT HAS TO BE RESTRICTED TO 75mm AS PER IS-CODES");
		    	fprintf(fq,"\n\n\nWidth of the footing is changed to %f", b);
                fprintf(fq,"\n\tIMMEDIATE SETTLEMENT OF SAND(in mm)= %f\n",S22);
				fprintf(fq,"\n\tAllowable Bearing Pressure(in kPa): %f",Qs1);
			}
			break;

		case 3:

        	Qs1=Q/A;
			Wel=t/nl;
			cen=Wel/2;
			for(i=1;i<=nl;i++)
			{
				if((cen+Df)<=Wl)
				{
					Po=(cen+Df)*y;
					P=Q/((b+cen)*(b+cen));
					FM[i]=log10((Po+P)/Po);
					TFM=TFM+FM[i];
					cen=cen+Wel;
	   			}
				else
				{
					Po=((y-9.81)*(Df+cen-Wl))+(y*Wl);
					P=Q/((b+cen)*(b+cen));
					FM[i]=log10((Po+P)/Po);
					TFM=TFM+FM[i];
					cen=cen+Wel;
				}
			}
			S1=1000*((Cc*Wel)/(1+e))*TFM;
			S2=1000*Qs1*b*((1-(m*m))/Es)*Is;
 		    S12=S1+S2;

	        cen1=(t1/2);
			if((cen1+Df+t)<=Wl)
			{
				Po=(cen1*y1)+((Df+t)*y);
				P=Q/((b+cen1+t)*(b+cen1+t));
				Val=log10((Po+P)/Po);
			}
			else if((Df+t)<=Wl)
			{
				Po=((y1-9.81)*(Df+t+cen1-Wl))+((y*(Df+t)))+(y1*(Wl-(Df+t)));
				P=Q/((b+cen1+t)*(b+cen1+t));
				Val=log10((Po+P)/Po);

			}

			else
			{
				Po=((y1-9.81)*cen1)+((y-9.81)*(Df+t-Wl))+(y*Wl);
				P=Q/((b+cen1+t)*(b+cen1+t));
				Val=log10((Po+P)/Po);
			}

		    Cr=400*N;
	    	Cs=((1.5*Cr)/Po);
			H1=t1;
		          S22=1000*(H1/Cs)*Val;

				  
			TFM=0;

			if(bstart==0)
			{
				fprintf(fq,"\n\nProposed Load on the Footing(in kN): %f",Q);
				fprintf(fq,"\n\tCONSOLIDATION SETTLEMENT OF CLAY(in mm)= %f",S1);
		    	fprintf(fq,"\n\tIMMEDIATE SETTLEMENT OF CLAY(in mm)= %f\n",S2);
			    fprintf(fq,"\n\tTOTAL SETTLEMENT OF CLAY(in mm)= %f\n",S12);
		    	fprintf(fq,"\n\tIMMEDIATE SETTLEMENT OF SAND(in mm)= %f",S22);
                fprintf(fq,"\n\tTOTAL SETTLEMENT OF SOIL (in mm)= %f\n",(S22+S12));
	            fprintf(fq,"\n\tAllowable Bearing Pressure(in kPa): %f",Qs1);
			}
			
			
			if((S12+S22)>=75)
			{
				bstart=1;
				if(bstart==1)
				b=b+0.05;
				goto bbegin;
				
			}

			if(bstart==1)
			{
				fprintf(fq,"\n\n**************************************************************");
				fprintf(fq,"\n\nAS THE SETTLEMENT HAS TO BE RESTRICTED TO 75mm AS PER IS-CODES");
			    fprintf(fq,"\n\nWidth of the footing is changed to %f", b);
				fprintf(fq,"\n\tCONSOLIDATION SETTLEMENT OF CLAY(in mm)= %f",S1);
		    	fprintf(fq,"\n\tIMMEDIATE SETTLEMENT OF CLAY(in mm)= %f\n",S2);
		    	fprintf(fq,"\n\tTOTAL SETTLEMENT OF CLAY(in mm)= %f\n",S12);
		    	fprintf(fq,"\n\tIMMEDIATE SETTLEMENT OF SAND(in mm)= %f",S22);
                fprintf(fq,"\n\tTOTAL SETTLEMENT OF SOIL (in mm)= %f\n",(S22+S12));
			
	        	 fprintf(fq,"\n\tAllowable Bearing Pressure(in kPa): %f",Qs1);	
			
			}
			break;

		case 4:

		    Qs1=Q/A;
	        cen1=(t1/2);

			if((cen1+Df)<=Wl)
			{
				Po=(cen1+Df)*y;
				P=Q/((b+cen1)*(b+cen1));
				Val=log10((Po+P)/Po);
			}
			else
			{
				Po=((y-9.81)*(Df+cen1-Wl))+(y*Wl);
				P=Q/((b+cen1)*(b+cen1));
				Val=log10((Po+P)/Po);
			}
			H1=t1;
			Cr=400*N;
			Cs=((1.5*Cr)/Po);
			         S22=1000*(H1/Cs)*Val;

		


            Wel=t/nl;
        	cen=Wel/2;

			for(i=1;i<=nl;i++)
			{
				if((cen+Df+t1)<=Wl)
				{
					Po=(cen*y1)+((Df+t1)*y);
					P=Q/((b+cen+t1)*(b+cen+t1));
					FM[i]=log10((Po+P)/Po);
					TFM=TFM+FM[i];
					cen=cen+Wel;
	   			}
				else if((Df+t1)<=Wl)
				{
					Po=((y1-9.81)*(Df+t1+cen-Wl))+((y*(Df+t1)))+(y1*(Wl-(Df+t1)));
					P=Q/((b+cen+t1)*(b+cen+t1));
					FM[i]=log10((Po+P)/Po);
					TFM=TFM+FM[i];
					cen=cen+Wel;
				}

				else
				{
					Po=((y1-9.81)*cen)+((y-9.81)*(Df+t1-Wl))+(y*Wl);
					P=Q/((b+cen+t1)*(b+cen+t1));
					FM[i]=log10((Po+P)/Po);
					TFM=TFM+FM[i];
					cen=cen+Wel;
				}
			}


	                S1=1000*((Cc*Wel)/(1+e))*TFM;
					
					S2=1000*Qs1*b*((1-(m*m))/Es)*Is;
					S12=S1+S2;

			TFM=0;
			if(bstart==0)
			{
			
				fprintf(fq,"\n\nProposed Load on the Footing(in kN): %f",Q);
				fprintf(fq,"\n\tIMMEDIATE SETTLEMENT OF SAND(in mm)= %f",S22);
	       		fprintf(fq,"\n\tCONSOLIDATION SETTLEMENT OF CLAY(in mm)= %f",S1);
 	            fprintf(fq,"\n\tIMMEDIATE SETTLEMENT OF CLAY(in mm)= %f\n",S2);
	            fprintf(fq,"\n\tTOTAL SETTLEMENT OF CLAY(in mm)= %f\n",S12);
	            fprintf(fq,"\n\tTOTAL SETTLEMENT OF SOIL (in mm)= %f\n",(S12+S22));
		    	fprintf(fq,"\n\tAllowable Bearing Pressure(in kPa): %f",Qs1);

			}
			
			
			if((S12+S22)>=75)
			{
				bstart=1;
				if(bstart==1)
				b=b+0.05;
				goto bbegin;
				
			}

			if(bstart==1)
			{
		    	fprintf(fq,"\n\n**************************************************************");
				fprintf(fq,"\n\nAS THE SETTLEMENT HAS TO BE RESTRICTED TO 75mm AS PER IS-CODES");
				fprintf(fq,"\n\nWidth of the footing is changed to %f", b);
		        fprintf(fq,"\n\tIMMEDIATE SETTLEMENT OF SAND(in mm)= %f",S22);
                fprintf(fq,"\n\tCONSOLIDATION SETTLEMENT OF CLAY(in mm)= %f",S1); 
	            fprintf(fq,"\n\tIMMEDIATE SETTLEMENT OF CLAY(in mm)= %f\n",S2);    
                fprintf(fq,"\n\tTOTAL SETTLEMENT OF CLAY(in mm)= %f\n",S12);
	            fprintf(fq,"\n\tTOTAL SETTLEMENT OF SOIL (in mm)= %f\n",(S12+S22));
			   	fprintf(fq,"\n\tAllowable Bearing Pressure(in kPa): %f",Qs1);
			}	
			break;

	default:

 	goto DL;
	break;

	}
	fclose(fp);
	fclose(fq);
	
}
