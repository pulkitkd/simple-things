/*    float K=400; //conductivity in W/m sq- K
	float h=15; //convective heat transfer coefficient in W/m-K
	float d=0.02; //diameter in m
	float L=0.3; //Length of fin in m
	float Tinf=25; //temperature of surrounding fluid in Degree Celsius
	float Tbase=150; //base temperature of fin in Degree Celsius
	int n=25; //no of nodes (1 more than no of divisions)
	float dx=L/(n-1); //grid size
	float A[n][n]; // coefficient matrix
	float X[n]; //variable matrix
	float B[n]; //constants

	for(int i=0;i<n;i++) //Initialize the matrix with all entries as zero
	for(int j=0;j<n;j++)
	A[i][j]=0;

	float P=3.14*d; //perimeter
	float Ar=3.14*(d*d)/4; //area of cross section

	float ae, aw, ap, b; // coefficients of discretized equation
	float Q=0; // flux at free end of fin

	ae=K/dx;
	aw=K/dx;
	ap=2*K/dx + h*P*dx/Ar ;
	b=h*P*Tinf*dx/Ar;

	//Boundary Conditions
	A[0][0]=1;
	B[0]=150;

	A[n-1][n-1]=K/dx + h*P*dx/Ar +h;
	A[n-1][n-2]=-aw;
	B[n-1]=b+h*Tinf;

	for (int i=1;i<=n-2;i++) //Entering values in coefficient matrix
	{
		A[i][i-1]=-aw;
		A[i][i]=ap;
		A[i][i+1]=-ae;
	}

	for(int i=0;i<n;i++) // display matrix
	{	for(int j=0;j<n;j++)
			cout<<A[i][j]<<"  ";
		cout<<endl;
	}

	float d1[n],d2[n],d3[n]; //arrays for diagonal elements, needed by TDMA solver

	for(int i=0;i<n;i++) //initialize the diagonal arrays
	{
		d1[i]=0;
		d2[i]=0;
		d3[i]=0;
	}

	for(int i=0;i<n-1;i++) //assign elements in the diagonal arrays
	d1[i+1]=A[i+1][i];

	for(int i=0;i<n;i++) //assign elements in the diagonal arrays
	d2[i]=A[i][i];

	for(int i=0;i<n-1;i++) //assign elements in the diagonal arrays
	d3[i]=A[i][i+1];

	for(int i=1;i<n-1;i++) //assign elements to array B; the array of constants
	B[i]=b;

	for(int i=0;i<n;i++) // display all the diagonals as column vectors, just a check
	{
		cout<<d1[i]<<" "<<d2[i]<<" "<<d3[i]<<" "<<B[i]<<endl;
	}

	solve(d1,d2,d3,B,n); // using the TDMA solver function 'solve'

	for (int i = 0; i < n; i++) // The B matrix resulting from 'solve' is the answer
		{
           cout << B[i] << endl;
        }

*/
