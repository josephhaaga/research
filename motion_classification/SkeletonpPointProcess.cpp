#include "stdafx.h"
#include "afxwin.h"
#include<string>
#include<vector>
#include<fstream>
#include<iostream>
#include<stdio.h>
using namespace std;
#pragma comment (lib, "winmm.lib" )

#include "SkeletonPointProcess.h"
//ofstream myfile;

//CFile myFile;
//CFileException e;
/*
//public:
	//double img_x, img_y;
//	double real_x, real_y, real_z;
	//double img_x_prev, img_y_prev;
CFile myFile;
	//  NP alter
//	std::vector<std::string>pathVector;
//	string coordinates[100];
	char buffer[50];
	// till here 
	
private:
	double absVel;
	bool initData;
	DWORD prevTime, elapTime;
*/

SkeletonPointProcess::SkeletonPointProcess()
	{
		
		this->Init();
		this->bPathVectorFill=false;
		//neil created 
		thetaLSprev = 0.0f;
		thetaRSprev = 0.0f;
		thetaLEprev = 0.0f;
		thetaREprev = 0.0f;

		thetaLSprev1 = 0.0f;
		thetaRSprev1 = 0.0f;
		thetaLEprev1 = 0.0f;
		thetaREprev1 = 0.0f;
		magposprev=0.0f;

		angvel_ls0=0.0f;
		angvel_rs0=0.0f;
		angvel_le0=0.0f;
		angvel_re0=0.0f;

		
		//neil end
		//myFile.Open (_T("example.txt"), CFile::modeCreate | CFile::modeWrite );
	}
SkeletonPointProcess::~SkeletonPointProcess()
{
	//myFile.Close();
}

std::vector<std::string>pathVector;
string coordinates[100];

	 void SkeletonPointProcess::Init()
	{
		this->initData=true;
		this->img_x = 0.0f;
		this->img_x_prev = 0.0f;
		this->img_y = 0.0f;
		this->img_y_prev = 0.0f;
		this->real_x = 0.0f;
		this->real_y = 0.0f;
		this->real_z = 0.0f;
		this->absVel = 0.0f;
		this->prevTime = timeGetTime();
		this->elapTime=1.0;
		this->timetrm=0.0f;
		this->printoutLS=0.0f;
		this->printoutRS=0.0f;
		this->printoutLE=0.0f;
		this->printoutRE=0.0f;
		this->printoutLS0=0.0f;
		this->printoutRS0=0.0f;
		this->printoutLE0=0.0f;
		this->printoutRE0=0.0f;
		this->motionunitnum=0.0f;
		this->weighttrm=0.0f;
		this->avgangleshouldhand=0.0f;
		this->avgangleneckelbow=0.0f;

		this->angxyplanehand=0.0f; //average
		this->angxyplanehandleft=0.0f;
		this->angxyplanehandright=0.0f;
		this->angxzplanehandleft=0.0f;
		this->angxzplanehandright=0.0f;
		this->angyzplanehandleft=0.0f;
		this->angyzplanehandright=0.0f;

		//this->
		
		//myfile.open ("example.txt");

		// NP alter
		//this->myfile.open("example.txt");
		for(int i= 0; i<100; i++)
			this->coordinates[i];// = new string[100]; // till here 
	}

	void SkeletonPointProcess:: SetVal(double imgx, double imgy, double realx, double realy, double realz)
	{
		this->img_x = imgx;
		this->img_y = imgy;
		this->real_x = realx;
		this->real_y = realy;
		this->real_z = realz;

		//



		
	}

	//void closeDataFile()
	//{

		//if(myfile.is_open());
		//{
		//	myfile.close();
		//}

		//till here 				
	//}// ADD VECTOR ADDITION FUNCTION

//	// NP alter:
//	void SkeletonPointProcess::addVector(double realx, double realy, double realz)
//	{
//		this->real_x = realx;
//		this->real_y = realy;
//		this->real_z = realz;
//
//		std::string x = to_string((long double)(realx));
//		std::string y = to_string((long double)(realy));
//		std::string z = to_string((long double)(realz));
//
//		//myfile.open ("example.txt");
//		//this->myFile.Open("c:\\Users_cfile_example.txt", CFile::modeCreate);
//		pathVector.push_back("(" + x + "," +  y + "," + z +")");
//
//		this->bPathVectorFill=true;
//
//	
////	myfile.open("example.txt", ios::out|ios::app);
////myfile <<"("+ x +"," + y + "," + z +")\n";
//		//myfile <<"("+ x +"," + y + "," + z +")\n";
//		
//	//	for(int i = 0 ; i<pathVector.size(); i++)
//	//	{
//
//	//	myfile <<"("+ x +"," + y + "," + z +")\n";
//	//	}
//	
//		
//		//closeDataFile();
//
//		//}
//	//	CString logstr;
//		//logstr.Format(_T("%s\n"),"("+ x + "," + y + ","+ z + ")");
//		//int iChars = 1+logstr.GetLength();
//	//	myFile.Write( (LPVOID)(LPCTSTR)logstr, iChars*sizeof(TCHAR) );
//		//myfile.open ("example.txt");
//		
//		
//
//	//	myfile<<"("+ x + y + z + ")";
//		//myFile.Write(" + x + y + z + ");
//		//myfile<<(x[0]);
//		
//	//	for(int i =0; i<pathVector.size();i++)
//		//{
//			
//			//myfile<<"("+ x +  y + z +")";
//		//}
//
//
//		
//		/*
//		CString logstr;
//		myFile.Open (_T("example.txt"), CFile::modeCreate | CFile::modeWrite );
//	
//		logstr.Format(_T("%s\n"),"("+ x + "," + y + ","+ z + ")");
//		int iChars = 1+logstr.GetLength();
//		myFile.Write( (LPVOID)(LPCTSTR)logstr, iChars*sizeof(TCHAR) );
//
//		myFile.Write<<"("+ x +"," + y + "," + z +")\n";
//		
//		myFile.Close();
//		*/
//		
//		
//
//
//	}
//	void SkeletonPointProcess:: copyVectorToArray()
//	{
//		if(pathVector.size()>0 && pathVector.size()<100)
//		{
//			for(int i=0; i<pathVector.size(); i++)
//			{
//				coordinates[i]= pathVector.at(i);
//			}
//		}
//		else if (pathVector.size()>=100)
//		{
//			int interval =pathVector.size()/100;
//			for (int i=0; i< 100 && interval*i<pathVector.size();i++)
//			{
//				coordinates[i]=pathVector.at(i*interval);
//
//			}
//		}
//	}
	//void closeDataFile()
	//	{
		//if(myfile.is_open()){
				//myfile.close();
		//}
	//}
		
	void SkeletonPointProcess::calcSpeed() //IMAGE BASED SPEED -> CHANGE LATER TO REAL SPEED
	{
		DWORD currTime = timeGetTime();
		double vx, vy;
		this->elapTime = currTime - this->prevTime;
	

		if(this->initData)
		{
			this->absVel = 0.0f;
			this->initData = false;
		}else{
			if(elapTime==0.0f)
			{
				this->absVel = 0.0f;
			}else{
				vx = (this->img_x - this->img_x_prev)/(double)elapTime;
				vy = (this->img_y - this->img_y_prev)/(double)elapTime;
				this->absVel = sqrt(vx*vx + vy*vy);
				newvel= sqrt(vx*vx + vy*vy);
			}
		}

		this->img_x_prev = this->img_x;
		this->img_y_prev = this->img_y;
		this->prevTime = currTime;
	}
	//Neil sample code try for time term
	void SkeletonPointProcess::calctimeterm(double thetaLS, double thetaRS, double thetaLE, double thetaRE)
	{
		//can do calctimetermandanglebetelbowhand
		int count=4;
		double time_ls=0;
		double time_rs=0;
		double time_le=0;
		double time_re=0;

		DWORD currTime = timeGetTime(); 
		this->elapTime = currTime - this->prevTime;
		//bottom to check for data existence
		if(this->initData)
		{
			this->absVel = 0.0f;
			this->initData = false;
		}else{
			if(elapTime==0.0f)
			{
				this->absVel = 0.0f;
			}else
			{//check for data ends//

				//note may have to change from degrees per second to radians per s. to do so, multiplye each weight_ls by 3.14/180
				
				time_ls=(thetaLS-thetaLSprev1)/elapTime;
				time_rs=(thetaRS-thetaRSprev1)/elapTime;
				time_le=(thetaLE-thetaLEprev1)/elapTime;
				time_re=(thetaRE-thetaREprev1)/elapTime;

				
				//***current angels equal previous angle as we move on to the next angle***
						thetaLSprev1=thetaLS;
						thetaRSprev1=thetaRS;
						thetaLEprev1=thetaLE;
						thetaREprev1=thetaRE;
				count=4;


				//do average of above two?
				timetrm=(180/PI)*(1000)*((time_ls+time_rs+time_le+time_re)/count); //this is the feature we are extracting; avg for all the joints
				
				//elbow
				//algorithm .... 
				

			}
		}
	}

	void SkeletonPointProcess::angreadingsfrom3dand2dvectorcalc(double thetaLS, double thetaRS, double thetaLE, double thetaRE, double thetaLS0, double thetaRS0, double thetaLE0, double thetaRE0)
	{


		DWORD currTime = timeGetTime(); 
		this->elapTime = currTime - this->prevTime;
		//bottom to check for data existence
		if(this->initData)
		{
			this->absVel = 0.0f;
			this->initData = false;
		}else{
			if(elapTime==0.0f)
			{
				this->absVel = 0.0f;
			}else
			{//check for data ends//


				printoutLS=thetaLS*(180/PI);
				printoutRS=thetaRS*(180/PI);
				printoutLE=thetaLE*(180/PI);
				printoutRE=thetaRE*(180/PI);
				printoutLS0=thetaLS0*(180/PI);
				printoutRS0=thetaRS0*(180/PI);
				printoutLE0=thetaLE0*(180/PI);
				printoutRE0=thetaRE0*(180/PI);
		


				//prof wanted non-averages for angles, but the below can be used for measures of sociability and outwardness
				//NOTE THESE TWO ARE BASED ON  THE 3D ANG VECTOR ALGORITHM
				avgangleneckelbow=(printoutLS+printoutRS)/2;
				avgangleshouldhand=(printoutLE + printoutRE)/2;

			}
		}
	}
	//	void SkeletonPointProcess::motionunits(double  thetaLS, double thetaRS, double thetaLE, double thetaRE)
	//{//for each joint calculate where absolute velocity is about zero
	////experiment with diff threshold values for abs velocity under which a pause can be signified and thus, one motion unit occurs
	//
	//	
	//	double d=sqrt(x^2+y^2+z^2)


	//	

	//	double thetaRSprev=0;
	//	double thetaLEprev=0;
	//	double thetaREprev=0;

	//	double angvel_ls0=0;
	//	double angvel_rs0=0;
	//	double angvel_le0=0;
	//	double angvel_re0=0;

	//	double angvel_ls1=0;
	//	double angvel_rs1=0;
	//	double angvel_le1=0;
	//	double angvel_re1=0;

	//	double thresholdvel=5.0f; //vary the threshold value for determining under which absolute velocity a motion unit occurs
	//	//*alternative approach: graph absVel vs time and approximate a threshold as a motion is carried out (with several motion units)*//

	//	DWORD currTime = timeGetTime(); 
	//	this->elapTime = currTime - this->prevTime;
	//	if(this->initData)
	//	{
	//		this->absVel = 0.0f;
	//		this->initData = false;
	//	}else{
	//		if(elapTime==0.0f)
	//		{

	//			this->absVel = 0.0f;
	//		}else
	//		{
	//			thetaLS=thetaLS*(180/3.14);
	//			thetaRS=thetaRS*(180/3.14);
	//			thetaLE=thetaLE*(180/3.14);
	//			thetaRE=thetaRE*(180/3.14);

	//		angvel_ls1=((thetaLS-thetaLSprev)/elapTime)*1000;
	//		angvel_rs1=((thetaRS-thetaRSprev)/elapTime)*1000;
	//		angvel_le1=((thetaLE-thetaLEprev)/elapTime)*1000;
	//		angvel_re1=((thetaRE-thetaREprev)/elapTime)*1000;

	//			//compare each angle

	//					thetaLSprev=thetaLS;
	//					thetaRSprev=thetaRS;
	//					thetaLEprev=thetaLE;
	//					thetaREprev=thetaRE;

	//			if (angvel_ls0>thresholdvel && angvel_ls1<thresholdvel)
	//			{
	//				motionunitnum=motionunitnum+1; //add a motion unit wherever the speed is 0; there is a pause in movement

	//			}
	//			
	//			if (angvel_rs0>thresholdvel && angvel_rs1<thresholdvel)
	//			{
	//				motionunitnum=motionunitnum+1; //add a motion unit wherever the speed is 0; there is a pause in movement

	//			}

	//			if (angvel_le0>thresholdvel && angvel_le1<thresholdvel)
	//			{
	//				motionunitnum=motionunitnum+1; //add a motion unit wherever the speed is 0; there is a pause in movement

	//			}

	//			if (angvel_re0>thresholdvel && angvel_re1<thresholdvel)
	//			{
	//				motionunitnum=motionunitnum+1; //add a motion unit wherever the speed is 0; there is a pause in movement

	//			}
	//			angvel_ls0=angvel_ls1;
	//			angvel_rs0=angvel_rs1;
	//			angvel_le0=angvel_le1;
	//			angvel_re0=angvel_re1;


	//			motionunitnum=angvel_le1;
	//		}
	//	}
	//}
	void SkeletonPointProcess::motionunitsandweight(double thetaLS, double thetaRS, double thetaLE, double thetaRE)
	{//for each joint calculate where absolute velocity is about zero
	//experiment with diff threshold values for abs velocity under which a pause can be signified and thus, one motion unit occurs
		
		//for weight
		double massLS=5000;
		double massRS=5000;
		double massLE=1000;
		double massRE=1000;
		double radLS=2;
		double radRS=2;
		double radLE=2;
		double radRE=2;
		//arbitrarily defined mass and radius values for joints considered

		double accelLS=0;
		double accelRS=0;
		double accelLE=0;
		double accelRE=0;
		double weight_ls=0;
		double weight_rs=0;
		double weight_le=0;
		double weight_re=0;
		double count=4;
		//for weight

		

		double angvel_ls1=0;
		double angvel_rs1=0;
		double angvel_le1=0;
		double angvel_re1=0;
		double diff=1.5; //vary the threshold value for determining under which absolute velocity a motion unit occurs
		//*alternative approach: graph absVel vs time and approximate a threshold as a motion is carried out (with several motion units)*//

		DWORD currTime = timeGetTime(); 
		this->elapTime = currTime - this->prevTime;
		if(this->initData)
		{
			this->absVel = 0.0f;
			this->initData = false;



		}else{
			if(elapTime==0.0f)
			{

				this->absVel = 0.0f;
			}else
			{
				thetaLS=thetaLS*(180/3.14);
				thetaRS=thetaRS*(180/3.14);
				thetaLE=thetaLE*(180/3.14);
				thetaRE=thetaRE*(180/3.14);

			angvel_ls1=((thetaLS-thetaLSprev)/elapTime)*1000;
			angvel_rs1=((thetaRS-thetaRSprev)/elapTime)*1000;
			angvel_le1=((thetaLE-thetaLEprev)/elapTime)*1000;
			angvel_re1=((thetaRE-thetaREprev)/elapTime)*1000;

				//compare each angle



				if ((angvel_le1-angvel_le0)>=diff && angvel_le0!=0)
				{
					motionunitnum=motionunitnum+1; //add a motion unit wherever the speed is 0; there is a pause in movement

				}
				
				if ((angvel_re1-angvel_re0)>=diff && angvel_re0!=0)
				{
					motionunitnum=motionunitnum+1; //add a motion unit wherever the speed is 0; there is a pause in movement

				}

				if ((angvel_ls1-angvel_ls0)>=diff && angvel_ls0!=0)
				{
					motionunitnum=motionunitnum+1; //add a motion unit wherever the speed is 0; there is a pause in movement

				}

				if ((angvel_rs1-angvel_rs0)>=diff && angvel_rs0!=0)
				{
					motionunitnum=motionunitnum+1; //add a motion unit wherever the speed is 0; there is a pause in movement

				}


			

				
				accelLS=(angvel_ls1-angvel_rs0)/elapTime;
				accelRS=(angvel_rs1-angvel_rs0)/elapTime;
				accelLE=(angvel_le1-angvel_le0)/elapTime;
				accelRE=(angvel_re1-angvel_re0)/elapTime;

				weight_ls=radLS*massLS*accelLS*sin((thetaLS-thetaLSprev));
				weight_rs=radRS*massRS*accelRS*sin((thetaRS-thetaRSprev));
				weight_le=radLE*massLE*accelLE*sin((thetaLE-thetaLEprev));
				weight_re=radRE*massRE*accelRE*sin((thetaRE-thetaREprev));
				weighttrm=(weight_ls+weight_rs+weight_le+weight_re)/count;

				angvel_ls0=angvel_ls1;
				angvel_rs0=angvel_rs1;
				angvel_le0=angvel_le1;
				angvel_re0=angvel_re1;




		
				
						thetaLSprev=thetaLS;
						thetaRSprev=thetaRS;
						thetaLEprev=thetaLE;
						thetaREprev=thetaRE;



			}
		}
	}

	//**************************THIS IS WRONG; IGNORE THIS BOTTOM PIECE OF CODE (THE ENTIRE FUNCTION), BUT DON'T COMMENT OUT YET NOR DELETE. 
	void SkeletonPointProcess::weightterm(double thetaLS, double thetaRS, double thetaLE, double thetaRE)
	{//for each joint, know: the radius of joint, its avg. mass, and use the input of acceleration and angle of that joint to determine torque which is radius*Force*sintheta=
		//radius*mass*acceleration*sin(theta)
	//

		double massLS=5000;
		double massRS=5000;
		double massLE=1000;
		double massRE=1000;
		double radLS=2;
		double radRS=2;
		double radLE=2;
		double radRE=2;
		//arbitrarily defined mass and radius values for joints considered - these may have to be adjusted to account for subject variations 

		//initializing initial joint angles to be zero 
		double thetaLSprev=0;
		double thetaRSprev=0;
		double thetaLEprev=0;
		double thetaREprev=0;
		//

		double accelLS=0;
		double accelRS=0;
		double accelLE=0;
		double accelRE=0;

		double angvel_ls0=0;
		double angvel_rs0=0;
		double angvel_le0=0;
		double angvel_re0=0;

		double angvel_ls1=0;
		double angvel_rs1=0;
		double angvel_le1=0;
		double angvel_re1=0;

		double weight_ls=0;
		double weight_rs=0;
		double weight_le=0;
		double weight_re=0;

		double count= 4; //let count represent the number of joints considered


		if(this->initData)
		{
			this->absVel = 0.0f;
			this->initData = false;
		}else{
			if(elapTime==0.0f)
			{
				this->absVel = 0.0f;
			}else
			{
				

				angvel_ls1=(thetaLS-thetaLSprev)/elapTime;
				angvel_rs1=(thetaRS-thetaRSprev)/elapTime;
				angvel_le1=(thetaLE-thetaLEprev)/elapTime;
				angvel_re1=(thetaRE-thetaREprev)/elapTime;


				
				accelLS=(angvel_ls1-angvel_rs0)/elapTime;
				accelRS=(angvel_rs1-angvel_rs0)/elapTime;
				accelLE=(angvel_le1-angvel_le0)/elapTime;
				accelRE=(angvel_re1-angvel_re0)/elapTime;

				weight_ls=radLS*massLS*accelLS*sin((thetaLS-thetaLSprev));
				weight_rs=radRS*massRS*accelRS*sin((thetaRS-thetaRSprev));
				weight_le=radLE*massLE*accelLE*sin((thetaLE-thetaLEprev));
				weight_re=radRE*massRE*accelRE*sin((thetaRE-thetaREprev));
								//***current angels equal previous angle as we move on to the next angle***
						thetaLSprev=thetaLS;
						thetaRSprev=thetaRS;
						thetaLEprev=thetaLE;
						thetaREprev=thetaRE;


				angvel_ls0=angvel_ls1;
				angvel_rs0=angvel_rs1;
				angvel_le0=angvel_le1;
				angvel_re0=angvel_re1;


				acttotalweighterm=(weight_ls+weight_rs+weight_le+weight_re)/count;
				}

			}
	}
	

	////Neil sample code ending

	double SkeletonPointProcess::GetAbsVel()
	{
		return this->absVel;
	}

	double SkeletonPointProcess::CalcAngleBtn3DVects(SkeletonPointProcess skp1, SkeletonPointProcess skp_origin,  SkeletonPointProcess skp2)
	{

		//neil:
		//do calculate linear velocity, linear acceleration, 
		double angle = 0.0f;

		double x1, y1, z1, x2, y2, z2;
		x1 = skp1.real_x - skp_origin.real_x;
		y1 = skp1.real_y - skp_origin.real_y;
		z1 = skp1.real_z - skp_origin.real_z;

		x2 = skp2.real_x - skp_origin.real_x;
		y2 = skp2.real_y - skp_origin.real_y;
		z2 = skp2.real_z - skp_origin.real_z;

		double len1 = sqrt( x1*x1 + y1*y1 + z1*z1 );
		double len2 = sqrt( x2*x2 + y2*y2 + z2*z2 );

		double dot_prod = ((x1*x2 + y1*y2 + z1*z2))/(len1*len2);
		double crs1 = y1*z2 - z1*y2;
		double crs2 = z1*x2 - x1*z2;
		double crs3 = x1*y2 - y1*x2;
		double cross_prod = (sqrt(crs1*crs1 + crs2*crs2 + crs3*crs3))/(len1*len2);

		angle = atan2(cross_prod, dot_prod);
		if(crs3<0)	angle = -angle;

		//double dot = x1*x2 + y1*y2 + z1*z2;

		//double val = dot / (len1 * len2);

		//if (val >= 1.0)
		//	angle = 0.0f;
		//else if (val <= -1.0)
		//	angle = 3.1415926/2.0f;
		//else
		//	angle = acos(val); // 0..PI




		return angle;
	}

	//my added function - neil 5/24/16 is below:
		double SkeletonPointProcess::CalcAngleBetweenplaneandHandparta(SkeletonPointProcess skp1_elb, SkeletonPointProcess skp2_hand, SkeletonPointProcess skp3_centershould, int planedet)
	{

		//neil:
		
		double angleret = 0.0f;
		
		double u1, u2, u3, A, B, C;
		u1 = skp2_hand.real_x - skp1_elb.real_x;
		u2 = skp2_hand.real_y - skp1_elb.real_y;
		u3 = skp2_hand.real_z - skp1_elb.real_z;
		
		
		if (planedet==1) //for xy plane
		{
		A=0;
		B=0;
		C=1;
		}

		else if (planedet==2) //for xz plane
		{
		A=0;
		B=1;
		C=0;
		}

		else if (planedet==3) //for yz plane
		{
		A=1;
		B=0;
		C=0;
		}
		double len1 = sqrt( A*A + B*B + C*C );
		double len2 = sqrt( u1*u1 + u2*u2 + u3*u3 );

		double sinw = (abs(A*u1 + B*u2 + C*u3))/(len1*len2);
	
	
		angleret = asin(sinw);
		
		//if dist or subtraction of hand coordinates minus 
		//if z of hand is less than z of center shoulder, than angle is negative  - need to implement this in bottom function
		//if ((skp2_hand.real_z - skp3_centershould.real_z)<=0)
		//	{angleret=-angleret;
		//	};

		return angleret;
	}

//do angles between yz and xz as well  
//				double SkeletonPointProcess::CalcAngleBetweenZYplaneandHandparta(SkeletonPointProcess skp1_elb, SkeletonPointProcess skp2_hand, SkeletonPointProcess skp3_centershould)
//	{
//
//		//neil:
//		
//		double angleret = 0.0f;
//		double sub=90.0f;
//
//		double u1, u2, u3, A, B, C;
//		u1 = skp2_hand.real_x - skp1_elb.real_x;
//		u2 = skp2_hand.real_y - skp1_elb.real_y;
//		u3 = skp2_hand.real_z - skp1_elb.real_z;
//
//		A=0;
//		B=0;
//		C=1;
//
//		double len1 = sqrt( A*A + B*B + C*C );
//		double len2 = sqrt( u1*u1 + u2*u2 + u3*u3 );
//
//		double sinw = (abs(A*u1 + B*u2 + C*u3))/(len1*len2);
//	
//	
//		angleret = asin(sinw);
//		
//		//if dist or subtraction of hand coordinates minus 
//		//if z of hand is less than z of center shoulder, than angle is negative  - need to implement this in bottom function
//		//if ((skp2_hand.real_z - skp3_centershould.real_z)<=0)
//		//	{angleret=-angleret;
//		//	};
//
//		return angleret;
//	}

		void SkeletonPointProcess::CalcAngleBetweenplaneandHandpartb(double thetaangleftxy, double thetaangrightxy,double thetaangleftxz, double thetaangrightxz,double thetaangleftyz, double thetaangrightyz, SkeletonPointProcess skp1_lh, SkeletonPointProcess skp2_rh, SkeletonPointProcess skp3_le, SkeletonPointProcess skp4_re)
	{

		//convert to degrees


		angxyplanehandleft=90-(thetaangleftxy*(180/PI));
		angxyplanehandright=90-(thetaangrightxy*(180/PI));
		angxzplanehandleft=90-(thetaangleftxz*(180/PI));
		angxzplanehandright=90-(thetaangrightxz*(180/PI));
		angyzplanehandleft=90-(thetaangleftyz*(180/PI));
		angyzplanehandright=90-(thetaangrightyz*(180/PI));
		
		///*important note: do remove bottom positive/negative sign determinant code IF too much noise and impractical and data is confused
		
		if ((skp1_lh.real_z - skp3_le.real_z)<0)
			{angxyplanehandleft=-angxyplanehandleft;
			};

		if ((skp2_rh.real_z - skp4_re.real_z)<0)
			{angxyplanehandright=-angxyplanehandright;
			};

		if ((skp1_lh.real_y - skp3_le.real_y)<0) //check orientation ie: is lh 
			{angxzplanehandleft=-angxzplanehandleft;
			};

		if ((skp2_rh.real_y - skp4_re.real_y)<0)
			{angxzplanehandright=-angxzplanehandright;
			};
		if ((skp1_lh.real_x - skp3_le.real_x)<0) //check orientation ie: is right hand with yz plane positive when should be negative 
			{angyzplanehandleft=-angyzplanehandleft;
			};

		if ((skp2_rh.real_x - skp4_re.real_x)<0)
			{angyzplanehandright=-angyzplanehandright;
			};


		//removal section ends
		angxyplanehand=(angxyplanehandright+angxyplanehandleft)/2; //not essential
		
	}

			////make a void function instance
//		double SkeletonPointProcess::calcbodyspeedacceler(SkeletonPointProcess skp_ex)
//	{
//				//can do calctimetermandanglebetelbowhand
//		
////speed and vel should be local valariables as we return them back 
//		
//		double magpos=0.0f;
//		double time_rs=0;
//
//		DWORD currTime = timeGetTime(); 
//		this->elapTime = currTime - this->prevTime;
//		//bottom to check for data existence
//		if(this->initData)
//		{
//			this->absVel = 0.0f;
//			this->initData = false;
//		}else{
//			if(elapTime==0.0f)
//			{
//				this->absVel = 0.0f;
//			}else
//			{//check for data ends//
//
//				magpos=sqrt((skp_ex.real_x*skp_ex.real_x) +(skp_ex.real_y*skp_ex.real_y)+(skp_ex.real_z*skp_ex.real_z));
//				vel_ex = (magpos-magposprev)/elapTime;
//				speed_ex=vel_ex;
//				if (vel_ex>=0)
//				{speed_ex=-speed_ex;
//				}
//
//				time_rs=(thetaRS-thetaRSprev1)/elapTime;
//			
//
//				
//				//***current angels equal previous angle as we move on to the next angle***
//						magposprev=magpos;
//						thetaRSprev1=thetaRS;
//
//				
//
//				
//		//neil:
//		//do calculate linear velocity, linear acceleration, 
//
//		return angle;
//	}

	////make a void function instance
	//	double SkeletonPointProcess::CalcAngleBtn3DVects(SkeletonPointProcess skp1, SkeletonPointProcess skp_origin,  SkeletonPointProcess skp2)
	//{

	//	//neil:
	//	//do calculate linear velocity, linear acceleration, 
	//	double angle = 0.0f;

	//	double x1, y1, z1, x2, y2, z2;
	//	x1 = skp1.real_x - skp_origin.real_x;
	//	y1 = skp1.real_y - skp_origin.real_y;
	//	z1 = skp1.real_z - skp_origin.real_z;

	//	x2 = skp2.real_x - skp_origin.real_x;
	//	y2 = skp2.real_y - skp_origin.real_y;
	//	z2 = skp2.real_z - skp_origin.real_z;

	//	double len1 = sqrt( x1*x1 + y1*y1 + z1*z1 );
	//	double len2 = sqrt( x2*x2 + y2*y2 + z2*z2 );

	//	double dot_prod = ((x1*x2 + y1*y2 + z1*z2))/(len1*len2);
	//	double crs1 = y1*z2 - z1*y2;
	//	double crs2 = z1*x2 - x1*z2;
	//	double crs3 = x1*y2 - y1*x2;
	//	double cross_prod = (sqrt(crs1*crs1 + crs2*crs2 + crs3*crs3))/(len1*len2);

	//	angle = atan2(cross_prod, dot_prod);
	//	if(crs3<0)	angle = -angle;

	//	//double dot = x1*x2 + y1*y2 + z1*z2;

	//	//double val = dot / (len1 * len2);

	//	//if (val >= 1.0)
	//	//	angle = 0.0f;
	//	//else if (val <= -1.0)
	//	//	angle = 3.1415926/2.0f;
	//	//else
	//	//	angle = acos(val); // 0..PI

	//	//neil added:
	//	



	//	return angle;
	//}

	double SkeletonPointProcess:: CalcAngleBtn2DVects(SkeletonPointProcess skp1, SkeletonPointProcess skp_origin,  SkeletonPointProcess skp2, int flag)// falg = 101, 110, 011
	{
		double angle = 0.0f;

		double x1, y1, x2, y2;
		//selecting among the plane in which to position the triangle
		if(flag==110){
			x1 = skp1.real_x - skp_origin.real_x;
			y1 = skp1.real_y - skp_origin.real_y;

			x2 = skp2.real_x - skp_origin.real_x;
			y2 = skp2.real_y - skp_origin.real_y;
		}else if(flag ==101){
			x1 = skp1.real_x - skp_origin.real_x;
			y1 = skp1.real_z - skp_origin.real_z;

			x2 = skp2.real_x - skp_origin.real_x;
			y2 = skp2.real_z - skp_origin.real_z;
		
		}else{//011
			x1 = skp1.real_y - skp_origin.real_y;
			y1 = skp1.real_z - skp_origin.real_z;
			
			x2 = skp2.real_y - skp_origin.real_y;
			y2 = skp2.real_z - skp_origin.real_z;
		}

		double len1 = sqrt( x1*x1 + y1*y1 ); //length 
		double len2 = sqrt( x2*x2 + y2*y2 );

		double dot_prod = (x1*x2 + y1*y2)/(len1*len2);
		//double val = dot / (len1 * len2);
		double crs1 = x1*y2 - y1*x2;
		double cross_prod = (abs(crs1))/(len2*len2);

		angle = atan2(cross_prod, dot_prod);
		if(crs1<0)	angle = -angle;

		//if (val >= 1.0)
		//	angle = 0.0f;
		//else if (val <= -1.0)
		//	angle = 3.1415926/2.0f;
		//else
		//	angle = acos(val); // 0..PI

		return angle;
	}
	
