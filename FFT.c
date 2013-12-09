float fft (float data[], usigned long number_of_complex_samples, int isign)
{
    n=number_of_complex_samples;
    
    for (i=0;i<n;i++) idx[i]=i;
    
    for (i=0;i<n/2;i+=2) {
        if (j > i) {
            idx[j] = i;
            idx[i] = j
            idx[j+1]=i+1;
            idx[i+1]=j+1;
            if((j/2)<(n/4)){
                idx[n-(i+2)]=n-(j+2);
                idx[n-(j+2)]=n-(i+2);
                idx[n-(i+2)+1]=n-(j+2)+1;
                idx[n-(j+2)+1]=n-(i+2)+1;
            }
        }
        m=n/2;
        while (m >= 2 && j >= m) {
            j -= m;
            m = m/2;
        }
        j += m;
    }
    
    idma_gather(v,data,idx,n,wordlen);
    
    mmax=2;
    while (mmax <= M*2)
    {
        istep = mmax<<1;
        theta=sinal*(2*pi/mmax);
        wtemp=sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi=sin(theta);
        wr=1.0;
        wi=0.0;
        for (m=1;m<mmax;m+=2) {
            for (i=m;i<=n;i+=istep) {
                j=i+mmax;
                tempr=wr*v[j-1]-wi*v[j];
                tempi=wr*v[j]+wi*v[j-1];
                v[j-1]=v[i-1]-tempr;
                v[j]=data[i]-tempi;
                v[i-1] += tempr;
                v[i] += tempi;
            }
            wr=(wtemp=wr)*wpr-wi*wpi+wr;
            wi=wi*wpr+wtemp*wpi+wi;
        }
        mmax=istep;
    }
    
    idma_transpose(v,v,m/M,M,ld1,ld2,wordlen);
    
    while (mmax < n)
    {
        istep = mmax<<1;
        theta=sinal*(2*pi/mmax);
        wtemp=sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi=sin(theta);
        wr=1.0;
        wi=0.0;
        for (m=1;m<mmax;m+=2) {
            for (i=m;i<=n;i+=istep) {
                j=i+mmax;
                tempr=wr*v[j-1]-wi*v[j];
                tempi=wr*v[j]+wi*v[j-1];
                v[j-1]=v[i-1]-tempr;
                v[j]=data[i]-tempi;
                v[i-1] += tempr;
                v[i] += tempi;
            }
            wr=(wtemp=wr)*wpr-wi*wpi+wr;
            wi=wi*wpr+wtemp*wpi+wi;
        }
        mmax=istep;
    }
}

float fft0 (float data[], unsigned long number_of_complex_samples, int isign)
{
    //variables for trigonometric recurrences
    unsigned long n,mmax,m,j,istep,i;
    double wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;

    //the complex array is real+complex so the array
    //as a size n = 2* number of complex samples
    // real part is the data[index] and
    //the complex part is the data[index+1]
    n=number_of_complex_samples * 2;
    
    //binary inversion (note that the indexes
    //start from 0 witch means that the
    //real part of the complex is on the even-indexes
    //and the complex part is on the odd-indexes
    j=0;
    for (i=0;i<n/2;i+=2) {
        if (j > i) {
            //swap the real part
            SWAP(data[j],data[i]);
            //swap the complex part
            SWAP(data[j+1],data[i+1]);
            // checks if the changes occurs in the first half
            // and use the mirrored effect on the second half
            if((j/2)<(n/4)){
                //swap the real part
                SWAP(data[(n-(i+2))],data[(n-(j+2))]);
                //swap the complex part
                SWAP(data[(n-(i+2))+1],data[(n-(j+2))+1]);
            }
        }
        m=n/2;
        while (m >= 2 && j >= m) {
            j -= m;
            m = m/2;
        }
        j += m;
    }
    //Danielson-Lanzcos routine
    mmax=2;
    //external loop
    while (n > mmax)
    {
        istep = mmax<<1;
        theta=sinal*(2*pi/mmax);
        wtemp=sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi=sin(theta);
        wr=1.0;
        wi=0.0;
        //internal loops
        for (m=1;m<mmax;m+=2) {
            for (i=m;i<=n;i+=istep) {
                j=i+mmax;
                tempr=wr*data[j-1]-wi*data[j];
                tempi=wr*data[j]+wi*data[j-1];
                data[j-1]=data[i-1]-tempr;
                data[j]=data[i]-tempi;
                data[i-1] += tempr;
                data[i] += tempi;
            }
            wr=(wtemp=wr)*wpr-wi*wpi+wr;
            wi=wi*wpr+wtemp*wpi+wi;
        }
        mmax=istep;
    }
}