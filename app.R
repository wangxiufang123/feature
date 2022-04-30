library(shiny)
library(dplyr)
###################functions###############################
cal_ncorr=function(x,y){
    n=length(x)
    dataf=data.frame(x=x,y=y)
    dataf_slice=arrange(dataf,x)
    dataf_slice[,2]=as.numeric(rank(dataf_slice[,2]))
    y1=dataf_slice[,2]
    s=sum(abs(dataf_slice[-1,2]-y1[-n]))
    ncorr=1-3*s/(n^2-1)
    return(ncorr)
}
cal_GMC=function(x,y,c){
    n=length(x)
    H=round(n/c)
    groupv=rep(0,H)
    dataf=data.frame(x=x,y=y)
    dataf_slice=arrange(dataf,x)
    for (h in 1:H){
        start=(h-1)*c+1
        end=h*c
        groupv[h]=var(dataf_slice[start:end,2])
    }
    gmc=1-mean(groupv)/var(y)
    return(gmc)
}
cal_hat_SIRS=function(x,y){
    n=length(y)
    x=scale(x, center = T,scale = T)
    y_copy=y
    w=rep(0,n)
    for (j in 1:n) {
        w[j]=(mean(x*as.numeric(y<y_copy[j])))^2
    }
    w_hat=mean(w)
    return(w_hat)
}
#DC_func
dcov=function(x,y){
    x_copy=x                    
    y_copy=y
    n=length(y)        
    S1=rep(0,n)
    S21=rep(0,n)
    S22=rep(0,n)
    S3=0
    for (i in 1:n) {      
        S1[i]=mean(abs(x-x_copy[i])*abs(y-y_copy[i]))
        S21[i]=mean(abs(x-x_copy[i]))
        S22[i]=mean(abs(y-y_copy[i]))
        for (j in 1:n) {
            S3=S3+sum(abs(x-x_copy[i])*abs(y-y_copy[j]))
        }
    }
    S11=mean(S1)
    S2=mean(S21)*mean(S22)
    S3=S3/(n^3)
    return(sqrt(S11+S2-2*S3))
}
cal_hat_dcorr=function(u,v){
    return(dcov(u,v)/sqrt(dcov(u,u))/sqrt(dcov(v,v)))
}
#MDC_func
cal_hat_MDC<-function(x,y){
    #construct matrix A and B
    d=length(x)
    A=matrix(0,ncol=d,nrow=d)
    Bm=matrix(0,ncol=d,nrow=d)
    for (i in 1:(d-1)){
        for (j in (i+1):d){
            A[i,j]=abs(x[i]-x[j])
            Bm[i,j]=(y[i]-y[j])^2
        }
    }
    A1=A+t(A)
    B1=(Bm+t(Bm))/2
    A2=matrix(0,ncol=d,nrow=d)
    B2=matrix(0,ncol=d,nrow=d)
    for (k in 1:(d-1)){
        for (l in k:d){
            A2[k,l]=A1[k,l]-rowSums(A1)[k]/(d-2)-colSums(A1)[l]/(d-2)+
                sum(A1)/(d-1)/(d-2)
            B2[k,l]=B1[k,l]-rowSums(B1)[k]/(d-2)-colSums(B1)[l]/(d-2)+
                sum(B1)/(d-1)/(d-2)
        }
    }
    A2_t=t(A2)
    diag(A2_t)=0
    A_tidle=A2+A2_t
    B2_t=t(B2)
    diag(B2_t)=0
    B_tidle=B2+B2_t
    #calculate test statistics Tn
    cn=(d-3)^4/(d-1)^4
    MDD_obs=2*sum(A2_t*B2_t)/d/(d-3)
    S2=2*sum(A2_t^2*B2_t^2)/cn/d/(d-1)
    Tn=sqrt(d*(d-1)/2)*MDD_obs/sqrt(S2)
    return(Tn)
}
#CD_func
cal_hat_CD=function(x,y){
    x_copy=x
    D_y=var(y)*(n-1)/n
    summ=0
    for (j in 1:length(x)){
        z=as.numeric(x<x_copy[j])
        summ=summ+(sum( (y-mean(y))*(z-mean(z)) ))^2
    }
    cd=summ/(length(x)^3)/D_y
    return(cd)
}

#############################################################
# Define UI for dataset viewer application
ui = shinyUI(pageWithSidebar(
    # Application title
    titlePanel("Feature Screening for High-dimensional Data"),
    
    sidebarPanel(
        helpText("Results of feature screening by different marginal dependence measures. "),
        
        #fileInput("dataset","Choose a csv dataset", accept = ".csv"),
        #textInput("response", "Input the response name:"),
        
        selectInput("method", label="Choose a marginal screening procedure:", 
                    choices = c( "SIS", "RRCS", "SIRS","DC-SIS","MDC-SIS", "CD-SIS","Sliced-GMC-SIS",'New-Corr-SIS'),
                    selected = "SIS"),
        
        numericInput("size", "Number of features to retain :", 10),
        numericInput("c", "Option: number of observations within each slice for Sliced-GMC-SIS:", 12,min = 2)
    ),
    
    
    mainPanel(
        h1("Screening results"),
        tableOutput("result"), 
        
    )
))

data=read.csv("./data0.csv")
coln=colnames(data)
y_index=which(coln=="1384040_at")
y=data[,y_index]
x=data[,-y_index]
p=dim(x)[2]


# Define server logic required to summarize and view the selected dataset
server = shinyServer(function(input, output) {
    #options(shiny.maxRequestSize=50*1024^2)
    
    output$caption <- renderText({
        input$caption
    })
    
    output$result <- renderTable({
        # data=input$dataset
        d=input$size
        if(input$method == 'SIS'){
            corr=as.numeric(abs(cor(x,y)))
            sel=which(rank(-corr)<=d)
            v2=rank(-corr)[sel]
            v3=corr[sel]
        }
        else if (input$method=='RRCS'){
            tau=as.numeric(abs(cor(x,y,method = "kendall")))
            sel=which(rank(-tau)<=d)
            v2=rank(-tau)[sel]
            v3=tau[sel]
        }
        else if(input$method=='SIRS'){
            SIRS=rep(0,p)
            for (l in 1:p) {
                SIRS[l]=cal_hat_SIRS(x[,l],y)
            }
            sel=which(rank(-SIRS)<=d)
            v2=rank(-SIRS)[sel]
            v3=SIRS[sel]
        }
        else if(input$method=='CD-SIS'){
            CD_bbs14=rep(0,p)
            for (l in 1:p) {
                CD_bbs14[l]=cal_hat_CD(x[,l],y)
            }
            sel=which(rank(-as.numeric(CD_bbs14))<=d)
            v2=rank(-as.numeric(CD_bbs14))[sel]
            v3=as.numeric(CD_bbs14)[sel]
        }
        else if(input$method=='DC-SIS'){
            DC_bbs14=rep(0,p)
            for (l in 1:p) {
                DC_bbs14[l]=cal_hat_dcorr(x[,l],y)
            }
            sel=which(rank(-as.numeric(DC_bbs14))<=d)
            v2=rank(-as.numeric(DC_bbs14))[sel]
            v3=as.numeric(DC_bbs14)[sel]
        }
        else if(input$method=='MDC-SIS'){
            MDC_bbs14=rep(0,p)
            for (l in 1:p) {
                MDC_bbs14[l]=cal_hat_MDC(x[,l],y)
            }
            sel=which(rank(-abs(MDC_bbs14))<=d)
            v2=rank(-abs(MDC_bbs14))[sel]
            v3=abs(MDC_bbs14)[sel]
        }
        else if(input$method=='Sliced-GMC-SIS'){
            c=input$c
            GMC_bbs14=rep(0,p)
            for (l in 1:p) {
                GMC_bbs14[l]=cal_GMC(x[,l],y,c)
            }
            GMC=abs(GMC_bbs14)
            sel=which(rank(-GMC)<=d)
            v2=rank(-GMC)[sel]
            v3=GMC[sel]
        }
        else if(input$method=='New-Corr-SIS'){
            ncorr_gene=rep(0,p)
            for (l in 1:p) {
                ncorr_gene[l]=cal_ncorr(x[,l],y)
            }
            ncorr=abs(ncorr_gene)
            sel=which(rank(-ncorr)<=d)
            v2=rank(-ncorr)[sel]
            v3=ncorr[sel]
        }
        
        v1=colnames(x)[sel]
        table0=as.data.frame(cbind(v2,v1,v3))
        table1=arrange(table0,desc(v3))
        colnames(table1)=c("rank","gene","corr")
        table1
        
    })
    
    
})

# Create Shiny app ----
shinyApp(ui = ui, server = server)


