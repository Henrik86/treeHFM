#include "SumProduct.h"



void sendMessageVToFHMT( factorGraph &fg, int j, int i){
    //std::vector<variableNode > v =fg.variableNodes;
    //std::vector<factorNode > f =fg.factorNodes;
    std::vector<double> msg_tmp = ones(fg.variableNodes.at(i).size);
    //%msg_tmp = ones(1,3);
    //%msg_tmp = 1;
    //cout<<"sendMessageVToFHMT i:"<<i<<" j: "<<j<<endl;
    double logS_tmp = 0.0;
    std::vector<int> diffInd = setdiff(fg.variableNodes.at(i).f_neighbors,j);
    if (j>i){
        
        if (fg.variableNodes.at(i).observed != 1){
            std::vector<int> diffInd = setdiff(fg.variableNodes.at(i).f_neighbors,j);
            for (int m=0; m<diffInd.size();m++){
                int k=diffInd.at(m);
                //for (k = setdiff(v(i).f_neighbors,j)){
                if (k <= i) {//% down f(k)-->x(i)
                    if(fg.factorNodes.at(k).msg_out_down.size()>0){
                        msg_tmp = multEl(msg_tmp , fg.factorNodes.at(k).msg_out_down.at(findMessage(fg.factorNodes.at(k).msg_out_down,i)).msg);
                        logS_tmp = logS_tmp + fg.factorNodes.at(k).acc_logS_msg_out_down.at(findMessage(fg.factorNodes.at(k).acc_logS_msg_out_down,i)).msg;
                    }
                }else {//% k > i up f(k)-->x(i)
                    if(fg.factorNodes.at(k).msg_out_up.size() >0){
                        msg_tmp = multEl(msg_tmp,fg.factorNodes.at(k).msg_out_up.at(findMessage(fg.factorNodes.at(k).msg_out_up,i)).msg);
                        logS_tmp = logS_tmp + fg.factorNodes.at(k).acc_logS_msg_out_up.at(findMessage(fg.factorNodes.at(k).acc_logS_msg_out_up,i)).msg;
                    }
                }
            }
            //
            double scale=sum(msg_tmp);
            //msg_tmp_norm=msg_tmp*(1/scale);
            //
            //fg.factorNodes.at(j).msg_in_down_scale.push_back(scale);
            msgD msg_in_down;
            msg_in_down.index=i;
            msg_in_down.msg=multEl(msg_tmp,(1/scale));
            //
            fg.factorNodes.at(j).msg_in_down.push_back(msg_in_down);
            fg.factorNodes.at(j).N_msg_in_down=fg.variableNodes.at(j).N_msg_in_down+1;
            msgAccLog mA;
            mA.index=i;
            mA.msg=logS_tmp + log(scale);
            fg.factorNodes.at(j).acc_logS_msg_in_down.push_back(mA); //% store the accumulated log scale
            //
            //cout<<"msg: "<<endl;
            //printVector(msg_in_down.msg);
            //cout<<"acc_logS_msg_in_down: "<<mA.msg<<endl;
            //cout<<"-----------------------"<<endl;
            //f(j).msg_in_down_norm=msg_tmp_norm;
        }else{
            msg_tmp = ones(fg.variableNodes.at(i).size);
            //msg_tmp = ones(1,3);
            //
            double scale=sum(msg_tmp);
            //msg_tmp_norm=msg_tmp*(1/scale);
            //
            //fg.factorNodes.at(j).msg_in_down_scale.push_back(scale);
            msgD msg_in_down;
            msg_in_down.index=i;
            msg_in_down.msg=multEl(msg_tmp,(1/scale));
            //
            fg.factorNodes.at(j).msg_in_down.push_back(msg_in_down);;
            fg.factorNodes.at(j).N_msg_in_down=fg.variableNodes.at(j).N_msg_in_down+1;
            msgAccLog mA;
            mA.index=i;
            mA.msg=log(scale);
            fg.factorNodes.at(j).acc_logS_msg_in_down.push_back(mA); //% store the accumulated log scale
            //cout<<"msg: "<<endl;
            //printVector(msg_in_down.msg);
            //cout<<"acc_logS_msg_in_down: "<<mA.msg<<endl;
            //cout<<"-----------------------"<<endl;
            //f(j).msg_in_down_norm = msg_tmp_norm;
        }
        
    }else{//i>j
        if (fg.variableNodes.at(i).observed != 1) {
            //for (k = setdiff(v(i).f_neighbors,j)){
            std::vector<int> diffInd = setdiff(fg.variableNodes.at(i).f_neighbors,j);
            for (int m=0; m<diffInd.size();m++){
                int k=diffInd.at(m);
                if(k <= i){ // down f(k)-->x(i)
                    //%for m=1:length(f(k).msg_out_down)
                    if(fg.factorNodes.at(k).msg_out_down.size()>0){
                        msg_tmp = multEl(msg_tmp,fg.factorNodes.at(k).msg_out_down.at(findMessage(fg.factorNodes.at(k).msg_out_down,i)).msg);
                        if(j==16 && i==17){
                            
                        }
                        double logValueMsg=0.0;
                        //if(fg.factorNodes.at(k).acc_logS_msg_out_down.size()>0){
                        logValueMsg=fg.factorNodes.at(k).acc_logS_msg_out_down.at(findMessage(fg.factorNodes.at(k).acc_logS_msg_out_down,i)).msg;
                        //}
                        logS_tmp = logS_tmp+logValueMsg;
                        
                    }
                    //%end
	               }else{ // k > i up f(k)-->x(i)
                       //%for m=1:length(f(k).msg_out_up)
                       if(fg.factorNodes.at(k).msg_out_up.size()>0){
                           msg_tmp = multEl(msg_tmp,fg.factorNodes.at(k).msg_out_up.at((findMessage(fg.factorNodes.at(k).msg_out_up,i))).msg);
                           if(j==16 && i==17){
                               
                           }
                           double logValueMsg=0.0;
                           //if(fg.factorNodes.at(k).acc_logS_msg_out_up.size()>0){
                           logValueMsg=fg.factorNodes.at(k).acc_logS_msg_out_up.at(findMessage(fg.factorNodes.at(k).acc_logS_msg_out_up,i)).msg;
                           //}
                           logS_tmp = logS_tmp + logValueMsg;
                       }
                       //%end
                   }
            }
            double scale=sum(msg_tmp);
            //%msg_tmp_norm=msg_tmp*(1/scale);
            //
            if(j==16 && i==17){
                
            }
            //fg.factorNodes.at(j).msg_in_up_scale.at(i) = scale;
            msgD msg_in_up;
            msg_in_up.index=i;
            msg_in_up.msg=multEl(msg_tmp,(1/scale));
            //
            fg.factorNodes.at(j).msg_in_up.push_back(msg_in_up);
            fg.factorNodes.at(j).N_msg_in_up=fg.variableNodes.at(j).N_msg_in_up+1;
            //%f(j).msg_in_up_norm=msg_tmp_norm;
            msgAccLog mA;
            mA.index=i;
            mA.msg=logS_tmp + log(scale);
            
            fg.factorNodes.at(j).acc_logS_msg_in_up.push_back(mA); //% store the accumulated log scale
            //cout<<"msg: "<<endl;
            //printVector(msg_in_up.msg);
            //cout<<"acc_logS_msg_in_up: "<<mA.msg<<endl;
            //cout<<"-----------------------"<<endl;
        }else{
            msg_tmp = ones(fg.variableNodes.at(i).size);
            //%msg_tmp = ones(1,3);
            //%
            double scale=sum(msg_tmp);
            //%msg_tmp_norm=msg_tmp*(1/scale);
            //%
            //fg.factorNodes.at(j).msg_in_up_scale.at(i) = scale;
            msgD msg_in_up;
            msg_in_up.index=i;
            msg_in_up.msg=multEl(msg_tmp,(1/scale));
            //
            fg.factorNodes.at(j).msg_in_up.push_back(msg_in_up);
            fg.factorNodes.at(j).N_msg_in_up=fg.variableNodes.at(j).N_msg_in_up+1;
            
            msgAccLog mA;
            mA.index=i;
            mA.msg=log(scale);
            fg.factorNodes.at(j).acc_logS_msg_in_up.push_back(mA);
            //cout<<"msg: "<<endl;
            //printVector(msg_in_up.msg);
            //cout<<"acc_logS_msg_in_up: "<<mA.msg<<endl;
            //cout<<"-----------------------"<<endl;
            //%f(j).msg_in_up_norm = msg_tmp_norm;
        }
    }
    fg.factorNodes.at(j).N_msg_in_sofar = fg.factorNodes.at(j).N_msg_in_sofar + 1;
    fg.variableNodes.at(i).N_msg_out_remain = fg.variableNodes.at(i).N_msg_out_remain - 1;
}
void sendMessageFToVHMT( factorGraph &fg, int j, int i,transitionMatrices tM){
    //UNTITLED3 Summary of this function goes here
    //   Send a message from factor node f(j) to variable node v(i)
    
    //std::vector<variableNode > v =fg.variableNodes;
    //std::vector<factorNode > f =fg.factorNodes;
    std::vector<double > msg_tmp = ones(fg.variableNodes.at(i).size);
    double logS_tmp = 0.0;
    //cout<<"sendMessageFToVHMT i:"<<i<<" j: "<<j<<endl;
    if(j>i){
        std::vector<double > msg_tmpS;
        if(fg.factorNodes.at(j).leaf != 1){
            int strain=0;
            std::vector<int> diffInd = setdiff(fg.factorNodes.at(j).v_neighbors,i);
            for (int m=0; m<diffInd.size();m++){
                int k=diffInd.at(m);
                if(k < j){//  down x(k)-->f(j)
                    msg_tmp =multEl(msg_tmp,fg.factorNodes.at(j).msg_in_down.at(findMessage(fg.factorNodes.at(j).msg_in_down,k)).msg);
                    logS_tmp = logS_tmp +fg.factorNodes.at(j).acc_logS_msg_in_down.at(findMessage(fg.factorNodes.at(j).acc_logS_msg_in_down,k)).msg;
                }else{// % k > i up x(k)-->f(j)
                    strain=k;
                    msg_tmp =multEl(msg_tmp,fg.factorNodes.at(j).msg_in_up.at(findMessage(fg.factorNodes.at(j).msg_in_up,k)).msg);
                    logS_tmp = logS_tmp + fg.factorNodes.at(j).acc_logS_msg_in_up.at(findMessage(fg.factorNodes.at(j).acc_logS_msg_in_up,k)).msg;
                    
                    
                }
            }
            
            std::vector<double > xM=msg_tmp;
            std::vector<double >xMN=msg_tmp;
            
            //replaces repmat because of speed
            //replaces repmat because of speed
            
            if(fg.factorNodes.at(j).division){
                strain=0;
                std::vector< std::vector<double > > msg_tmp_strainUp;
                std::vector<int> diffInd = setdiff(fg.factorNodes.at(j).v_neighbors,i);
                for (int m=0; m<diffInd.size();m++){
                    int k=diffInd.at(m);
                    msg_tmp_strainUp.push_back(fg.factorNodes.at(j).msg_in_up.at(findMessage(fg.factorNodes.at(j).msg_in_up,k)).msg) ;
                }
                int variableIndex=0; // i do not need it as everything is summed
                //if(strain>i)
                // variableIndex=2;
                //end
                //UP
                std::vector<int>v_neighbors=fg.factorNodes.at(j).v_neighbors;
                int minIndex=min(fg.factorNodes.at(j).v_neighbors).index;
                v_neighbors.erase(v_neighbors.begin() + minIndex);
                int minN=min(v_neighbors).value;
                std::vector<double >xEU=fg.factorNodes.at(j).probO.at(minN); //'index on + strain is smaller
                
                int maxN=max(v_neighbors).value;
                std::vector<double >xED=fg.factorNodes.at(j).probO.at(maxN); //'index on - strain is larger
                //+ strain
                std::vector<double >oMU;
                for(int c=0;c<fg.variableNodes.at(i).size;c++){
                    for (int l=0;l<fg.variableNodes.at(i).size;l++){
                        oMU.push_back(xEU.at(c));
                    }
                }
                //- strain
                std::vector<double >oMD;
                for(int c=0;c<fg.variableNodes.at(i).size;c++){
                    for (int l=0;l<fg.variableNodes.at(i).size;l++){
                        oMD.push_back(xED.at(c));
                    }
                }
                oMD=arrangeItems(oMD,fg.factorNodes.at(j).indicesXD); //rearrange the entries for X2
                //UP
                std::vector<double> xMU;
                for(int c=0;c<fg.variableNodes.at(i).size;c++){
                    for (int l=0;l<fg.variableNodes.at(i).size;l++){
                        xMU.push_back(msg_tmp_strainUp.at(0).at(c));
                    }
                }
                
                //DOWN
                std::vector<double> xMD;
                for (int c=0;c<fg.variableNodes.at(i).size;c++){
                    for (int l=0;l<fg.variableNodes.at(i).size;l++){
                        xMD.push_back(msg_tmp_strainUp.at(1).at(c));
                    }
                }
                xMD=arrangeItems(xMD,fg.factorNodes.at(j).indicesXD); //rearrange the entries for X2
                //
                std::vector<std::vector<double > > xMUrep=repmat(xMU,tM.probDiv.at(variableIndex).size(),1);
                std::vector<std::vector<double > > xMDrep=repmat(xMD,tM.probDiv.at(variableIndex).size(),1);
                std::vector<std::vector<double > > oMUrep= repmat(oMU,tM.probDiv.at(variableIndex).size(),1);
                std::vector<std::vector<double > > oMDrep= repmat(oMD,tM.probDiv.at(variableIndex).size(),1);
                //
                std::vector<std::vector<double > >msg_tmpM;
                msg_tmpM=multEl(xMUrep,oMUrep);
                if(j==21 && i==20){
                    
                }
                
                msg_tmpM=multEl(msg_tmpM,tM.probDiv.at(variableIndex));
                if(j==21 && i==20){
                    
                }
                msg_tmpM=multEl(msg_tmpM,oMDrep);
                if(j==21 && i==20){
                }
                msg_tmpM=multEl(msg_tmpM,xMDrep);
                if(j==21 && i==20){
                }
                msg_tmpS=sum(msg_tmpM,1);
                
                
            }else{
                
                std::vector<double >xE=fg.factorNodes.at(j).probO.at(max(fg.factorNodes.at(j).v_neighbors).value);//';
                
                ///
                std::vector<std::vector<double > >xErep=repmat(xE,tM.probSeq.at(0).size(),1);
                std::vector<std::vector<double > >xMrep=repmat(xM,tM.probSeq.at(0).size(),1);
                
                std::vector<std::vector<double > >msg_tmpM=multEl(xMrep,xErep);
                
                msg_tmpM=multEl(msg_tmpM,tM.probSeq.at(0));
                msg_tmpS=sum(msg_tmpM,1);
                
            }
            
            
        }
        // transpose message
        //
        double scale=sum(msg_tmpS);
        //
        //
        
        msgD msg_out_up;
        msg_out_up.index=i;
        msg_out_up.msg=multEl(msg_tmpS,(1/scale));
        //cout<<"out up:"<<endl;
        //printVector(msg_out_up.msg);
        fg.factorNodes.at(j).msg_out_up.push_back(msg_out_up);
        fg.factorNodes.at(j).msg_out_up_scale.push_back(scale);
        msgAccLog mA;
        mA.index=i;
        mA.msg=logS_tmp + log(scale);
        
        fg.factorNodes.at(j).acc_logS_msg_out_up.push_back(mA);
        //cout<<"msg: "<<endl;
        //printVector(msg_out_up.msg);
        //cout<<"acc_logS_msg_out_up: "<<mA.msg<<endl;
    }else{
        //j;
        int obs=0;
        std::vector<double >msg_tmpM_summed;
        if (fg.factorNodes.at(j).leaf != 1){
            std::vector<int> diffInd = setdiff(fg.factorNodes.at(j).v_neighbors,i);
            double logS_tmp = 0.0;
            for (int m=0; m<diffInd.size();m++){
                int k=diffInd.at(m);
                if(k < j){ //% down x(k)-->f(j)
                    msg_tmp = multEl(msg_tmp,fg.factorNodes.at(j).msg_in_down.at(findMessage(fg.factorNodes.at(j).msg_in_down,k)).msg);
                    logS_tmp = logS_tmp + fg.factorNodes.at(j).acc_logS_msg_in_down.at(findMessage(fg.factorNodes.at(j).acc_logS_msg_in_down,k)).msg;
                }else{ //% k > i up x(k)-->f(j)
                    msg_tmp = multEl(msg_tmp,fg.factorNodes.at(j).msg_in_up.at(findMessage(fg.factorNodes.at(j).msg_in_up,k)).msg);
                    logS_tmp = logS_tmp + fg.factorNodes.at(j).acc_logS_msg_in_up.at(findMessage(fg.factorNodes.at(j).acc_logS_msg_in_up,k)).msg;
                }
            }
            if(fg.factorNodes.at(j).division){ //check if division =msg_tmp*(1/scale) true and summ states
                //cout<<"division"<<endl;
                std::vector<double >msg_tmp_strainUp = ones(fg.variableNodes.at(i).size);
                std::vector<double >msg_tmp_strainDown= ones(fg.variableNodes.at(i).size);
                int strain=0;
                std::vector<int> diffInd = setdiff(fg.factorNodes.at(j).v_neighbors,i);
                for (int m=0; m<diffInd.size();m++){
                    int k=diffInd.at(m);
                    if(k < j){ // down x(k)-->f(j)
                        msg_tmp_strainDown = multEl(msg_tmp_strainDown,fg.factorNodes.at(j).msg_in_down.at(findMessage(fg.factorNodes.at(j).msg_in_down,k)).msg);
                    }else{ // k > i up x(k)-->f(j)
                        strain=k;
                        msg_tmp_strainUp = multEl(msg_tmp_strainUp,fg.factorNodes.at(j).msg_in_up.at(findMessage(fg.factorNodes.at(j).msg_in_up,k)).msg);
                    }
                }
                int variableIndex=0;
                if(strain<i){
                    variableIndex=1;
                }
                std::vector<int>v_neighbors=fg.factorNodes.at(j).v_neighbors;
                Min minNN = min(fg.factorNodes.at(j).v_neighbors);
                int minIndex=minNN.index;
                v_neighbors.erase(v_neighbors.begin() + minIndex);
                Min minV =min(v_neighbors);
                int minN=minV.value;
                
                std::vector<double >xEU=fg.factorNodes.at(j).probO.at(minN); //'index on + strain is smaller
                Max maxV=max(v_neighbors);
                int maxN=maxV.value;
                std::vector<double >xED=fg.factorNodes.at(j).probO.at(maxN); //'index on - strain is larger
                //+ strain
                std::vector<double >oMU;
                for(int c=0;c<fg.variableNodes.at(i).size;c++){
                    for (int l=0;l<fg.variableNodes.at(i).size;l++){
                        oMU.push_back(xEU.at(c));
                    }
                }
                //- strain
                std::vector<double > oMD;
                for(int c=0;c<fg.variableNodes.at(i).size;c++){
                    for (int l=0;l<fg.variableNodes.at(i).size;l++){
                        oMD.push_back(xED.at(c));
                    }
                }
                
                oMD=arrangeItems(oMD,fg.factorNodes.at(j).indicesXD); //rearrange the entries for X2
                std::vector<std::vector<double > > xMDown=repmat(msg_tmp_strainDown,tM.probDiv.at(variableIndex).at(0).size(),1);
                std::vector<std::vector<double > > xMDownT=transpose(xMDown);
                //UP
                std::vector<double > xMUp;
                for(int c=0;c<tM.probDiv.at(variableIndex).size();c++){
                    for (int l=0;l<msg_tmp_strainUp.size();l++){
                        xMUp.push_back(msg_tmp_strainUp.at(l));
                    }
                }
                std::vector<std::vector<double > > xMUprep=repmat(xMUp,tM.probDiv.at(variableIndex).size(),1);
                //
                std::vector<std::vector<double > > oMUrep=repmat(oMU,tM.probDiv.at(0).size(),1);
                
                std::vector<std::vector<double > > oMDrep=repmat(oMD,tM.probDiv.at(0).size(),1);
                std::vector<std::vector<double > > transmatObs=multEl(oMUrep,tM.probDiv.at(0));
                //
                transmatObs=multEl(transmatObs,oMDrep);
                if(variableIndex==1){
                    transmatObs=arrangeItems(transmatObs,fg.factorNodes.at(j).indicesXD);
                }
                
                std::vector<std::vector<double > > msg_tmpM;
                msg_tmpM=multEl(xMDownT,transmatObs);
                
                msg_tmpM=multEl(msg_tmpM,xMUprep);
                
                //
                std::vector<double > msg_tmp = sum(msg_tmpM,2);
                
                std::vector<std::vector<int> > indicesToSum=fg.factorNodes.at(j).summationIndices;
                //
                //
                //
                for(int n=0; n<indicesToSum.size();n++){
                    msg_tmpM_summed.push_back(sum(arrangeItems(msg_tmp,indicesToSum.at(n))));
                }
            }else{
                std::vector<double > xM=msg_tmp;
                std::vector<std::vector<double > > xMrep=repmat(xM,tM.probDiv.at(0).size(),1);
                std::vector<std::vector<double > > xMrepT=transpose(xMrep);
                std::vector<double> xE=fg.factorNodes.at(j).probO.at(i); //' format [p(o|r) p(o|g) ...]
                std::vector<std::vector<double > >xErep=repmat(xE, tM.probDiv.at(0).size(),1);
                
                obs=0;
                std::vector<std::vector<double > > msg_tmpM;
                msg_tmpM=multEl(xMrepT,xErep);
                msg_tmpM=multEl(msg_tmpM,tM.probSeq.at(0));
                msg_tmpM_summed=sum(msg_tmpM,2);
            }
        }else{
            std::vector<double> xE=fg.factorNodes.at(j).probO.at(i);
            msg_tmp = multEl(fg.factorNodes.at(j).msg_in_down.at(findMessage(fg.factorNodes.at(j).msg_in_down,1)).msg,tM.probStart.at(0).at(0));
            msg_tmpM_summed = multEl(msg_tmp,xE);
        }
        //
        
        double scale=sum(msg_tmpM_summed);
        if(obs==1){
            scale=1;
        }
        //
        //
        msgD msg_out_down;
        msg_out_down.index=i;
        msg_out_down.msg=multEl(msg_tmpM_summed,(1/scale));
        fg.factorNodes.at(j).msg_out_down.push_back(msg_out_down);
        
        fg.factorNodes.at(j).msg_out_down_scale.push_back(scale);
        msgAccLog mA;
        mA.index=i;
        mA.msg=logS_tmp + log(scale);
        fg.factorNodes.at(j).acc_logS_msg_out_down.push_back(mA);
        //cout<<"msg: "<<endl;
        //printVector(msg_out_down.msg);
        //cout<<"acc_logS_msg_out_down: "<<mA.msg<<endl;
    }
    fg.variableNodes.at(i).N_msg_in_sofar = fg.variableNodes.at(i).N_msg_in_sofar + 1;
    fg.factorNodes.at(j).N_msg_out_remain = fg.factorNodes.at(j).N_msg_out_remain - 1;
}



factorGraph sumProductAlgorithm(factorGraph fgS, std::vector<int> endNodes, std::vector<int> startNodes, transitionMatrices tM){
    factorGraph fg=fgS;
    
    int done = 0;
    //% logS = 0; % log of the scale factor S, where S is the normalization factor of each message
    int counter=0;
    int counterM=0;
    //%%%message passing upwards
    string msg;
    while (endNodes.size()>0){
        int actNode=endNodes.at(0);
        //%send messages from the end node until it is not possible any
        //%more
        int stopSending=0;
        if(fg.variableNodes.at(actNode).f_children.size()==0){
            msg="v-->f"; //%ending variable node of a strain
        }else{
            msg="f-->v"; //%division factor node
            if (fg.factorNodes.at(actNode).N_msg_in_sofar >= fg.factorNodes.at(actNode).N_v_neighbors - 1){ //%the node was already visited
                stopSending=1;
                endNodes.erase(endNodes.begin());
            }
        }
        
        while(stopSending==0){
            
            if(msg=="f-->v"){
                std::vector<int> v_parent=fg.factorNodes.at(actNode).v_parents;
                
                if(v_parent.size()==0){
                    if(intersect(endNodes,actNode)==0){
                        endNodes.push_back(actNode);//%%add the actual node to the stack
                    }
                    endNodes.erase(endNodes.begin());//%drop the endNode from the stack
                    stopSending=1;
                    break;
                }else{
                    if (fg.factorNodes.at(actNode).N_msg_in_sofar >= fg.factorNodes.at(actNode).N_v_neighbors - 1){ //if there have not enough messages arrived, stop go to the next end node
                        sendMessageFToVHMT(fg,actNode,v_parent.at(0),tM); //only one v_parent possible
                    }else{
                        if(intersect(endNodes,actNode)==0){
                            endNodes.push_back(actNode);//%%add the actual node to the stack
                        }
                        endNodes.erase(endNodes.begin());//%drop the endNode from the stack
                        stopSending=1;
                        break;
                    }
                    actNode=v_parent.at(0); //only one v_parent possible
                    msg="v-->f";
                }
            }else if(msg=="v-->f"){
                if(fg.variableNodes.at(actNode).f_neighbors.size()==0){
                    if(intersect(endNodes,actNode)==0 ){
                        endNodes.push_back(actNode);//%%add the actual node to the stack
                    }
                    endNodes.erase(endNodes.begin());//%drop the endNode from the stack
                    stopSending=1;
                    break;
                }else{
                    int f_parent=min(fg.variableNodes.at(actNode).f_neighbors).value;
                    
                    sendMessageVToFHMT(fg,f_parent,actNode);
                    actNode=f_parent;
                    msg="f-->v";
                }
            }
        }
    }
    //%%%message passing downwards
    //cout<<"----------DOWN-------------"<<endl;
    while(startNodes.size()>0){
        
        int actNode=startNodes.at(0);
        //%send messages from the end node until it is not possible any
        //%more
        int stopSending=0;
        msg="f-->v"; //%start node is always a factor node
        while (stopSending==0){
            if(msg=="f-->v"){
                std::vector<int> v_child=fg.factorNodes.at(actNode).v_children;
                if(v_child.size()>0){ //%%if there have not enough messages arrived, stop go to the next end node
                    int v_childS;
                    if(v_child.size()>1){//%%division
                        if(intersect(startNodes,actNode)==0){
                            startNodes.push_back(actNode); //add the other strain to the startNodes
                            v_childS=max(v_child).value;
                        }else{ //%+ strain was already visited, start with - strain
                            v_childS=min(v_child).value;
                        }
                    }else{
                        v_childS=v_child.at(0);
                    }
                    
                    sendMessageFToVHMT(fg,actNode,v_childS,tM);
                    actNode=v_childS;
                    msg="v-->f";
                }else{
                    startNodes.erase(startNodes.begin());//%drop the startNode from the stack
                    stopSending=1;
                    break;
                }
                
            }else if(msg=="v-->f"){
                if(fg.variableNodes.at(actNode).f_children.size()>0){
                    int f_child=max(fg.variableNodes.at(actNode).f_children).value;
                    sendMessageVToFHMT(fg,f_child,actNode);
                    actNode=f_child;
                    msg="f-->v";
                }else{
                    startNodes.erase(startNodes.begin());//drop the startNode from the stack
                    stopSending=1;
                    break;
                }
            }
            
        }
    }
    //cout<<"SM ended"<<endl;
    return fg;
}



