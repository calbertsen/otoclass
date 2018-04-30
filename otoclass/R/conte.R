## From Morphometrics with R
## Modified to start at right most non-zero pixel
Conte<-function(imagematrix){
    indx <- which(imagematrix > 0,arr.ind=TRUE)
    x <- indx[which.max(indx[,1])[1],]
    I<-imagematrix
    ## x<-rev(x)
    ## x[1]<-dim(I)[1]-x[1]
    while (abs(I[x[1],x[2]]-I[x[1],(x[2]-1)])<0.1){x[2]<-x[2]-1}
    a<-1
    M<-matrix(c(0,-1,-1,-1,0,1,1,1,1,1,0,-1,-1,-1,0,1),
              2,8,byrow=T)
    M<-cbind(M[,8],M,M[,1])
    X<-0; Y<-0;
    x1<-x[1]; x2<-x[2]
    SS<-NA; S<-6
    while ((any(c(X[a],Y[a])!=c(x1,x2) ) | length(X)<3))
    {if (abs(I[x[1]+M[1,S+1],x[2]+M[2,S+1]]-I[x[1],x[2]])<0.1)
     {a<-a+1;X[a]<-x[1];Y[a]<-x[2];x<-x+M[,S+1]
         SS[a]<-S+1; S<-(S+7)%%8}
     else if (abs(I[x[1]+M[1,S+2],x[2]+M[2,S+2]]
                  -I[x[1],x[2]])<0.1)
        {a<-a+1;X[a]<-x[1];Y[a]<-x[2];x<-x+M[,S+2]
            SS[a]<-S+2; S<-(S+7)%%8}
     else if (abs(I[x[1]+M[1,(S+3)],x[2]+M[2,(S+3)]]
                  -I[x[1],x[2]])<0.1)
        {a<-a+1;X[a]<-x[1];Y[a]<-x[2];x<-x+M[,(S+3)]
            SS[a]<-S+3; S<-(S+7)%%8}
     else S<-(S+1)%%8}
    list(X=(Y[-1]), Y=((dim(I)[1]-X))[-1])
}
