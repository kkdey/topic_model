
##############  Clustering matrix by Factorization ########################

####  The goal of this script is to find the approximation C=WH, we shall

####  be given the C matrix or the counts matrix and we have to solve for 

####  W and H matrices.
N=dim(omega_true)[1]; G=dim(alpha_true)[2]; 
H_true=alpha_true+1000;
read_counts=omega_true%*%H_true +0*matrix(rnorm(N*G,0,1),nrow=N);

counts=read_counts;  

scale=3; K=4;

omega0=matrix(rdirichlet(N.pairs*T,c(scale/K,scale/K,scale/K,scale/K)), nrow=N.pairs*T);

iteration=1;
MaxIter=500
diff2=1000; diff1=100;

while(iteration < MaxIter)
{

	####   Estimation of the matrix H

	svd_omega=svd(omega0);
	temp1=t(svd_omega$v)%*%diag(1/svd_omega$d^2,dim(omega0)[2])%*%svd_omega$v;
	temp2=t(omega0)%*%counts;
	temp1=solve(t(omega0)%*%omega0);
	H = temp1%*%temp2;

	###  Estimation of the matrix W (or omega) 

	omega=matrix(0,dim(counts)[1],K);
	for(n in 1:dim(counts)[1])
	{
		omega_vec=omega0[n,];
		counts_vec=counts[n,];
		res=optim(reverse_transform(omega_vec), function(v) loglik_norm(transform(v),counts_vec,t(H)) );
		omega[n,]=transform(res$par);
	}
	diff2=diff1;
	diff1=fnorm(counts,omega%*%H);
	cat("The difference is",diff1,"\n");
	omega0=omega;
	iteration=iteration+1;
}

	docweights=omega0;
	library(permute);
	library("BioPhysConnectoR");
	perm_set=rbind(1:K,allPerms(1:K));
	diff=array(0,dim(perm_set)[1]);
	for (p in 1:dim(perm_set)[1])
	{
		temp=docweights[,perm_set[p,]];
		diff[p]=fnorm(temp,omega_true);
	}

	p_star=which(diff==min(diff));
	docweights=docweights[,perm_set[p_star,]];


	windows()
	par(mar=c(8,5.25,2.5,2.5))

	# - get rid of space space between leftmost bar and y axis
	par(xaxs="i")

	k=K
	# Make plot 
	# - plot data
	barplot(t(docweights),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)

	windows()
	barplot(t(omega_true),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)


	############  usual topic model  fit  ##########################


	Topic_Clus=topics(counts,K,kill=0,tol=0.01);
	docweights_topics=Topic_Clus$omega;
	library(permute);
	library("BioPhysConnectoR");
	perm_set=rbind(1:K,allPerms(1:K));
	diff=array(0,dim(perm_set)[1]);
	for (p in 1:dim(perm_set)[1])
	{
		temp=docweights_topics[,perm_set[p,]];
		diff[p]=fnorm(temp,omega_true);
	}

	p_star=which(diff==min(diff));
	docweights_topics=docweights_topics[,perm_set[p_star,]];



	windows()
	par(mar=c(8,5.25,2.5,2.5))

	# - get rid of space space between leftmost bar and y axis
	par(xaxs="i")

	barplot(t(docweights_topics),col=2:(k+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)


