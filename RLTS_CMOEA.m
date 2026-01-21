    classdef RLTS_CMOEA < ALGORITHM
    % <multi> <real/integer/label/binary/permutation> <constrained>
    
    
        methods
            function main(Algorithm,Problem)
                Population1 = Problem.Initialization();
                Population2 = Problem.Initialization(); 
                Population3 = Problem.Initialization(); 
                Population = [Population1,Population2,Population3];
               
                [RW,WVs] = UniformPoint(Problem.N,Problem.M);
                T = ceil(WVs/10);
                nr = ceil(WVs/100);
                W = zeros(size(RW));
                for i = 1:WVs
                    W(i, :) = RW(WVs - i + 1, :);
                end
                
                B = pdist2(W,W);
                [~,B] = sort(B,2);
                B = B(:,1:T);
    
                Q = zeros(2,size(WVs,1));
                M = Q;
                flag2 = 0;
                alp = 0.1;
                gam = 0.95;
                greedy = 0.1;
                epsilon = 0.1;
                test_Q2 = randi(size(WVs,1), 1, Problem.N);
                angle=acos(1-pdist2(W,W,'cosine'));
                temp_angle=angle;
                temp_angle(logical(eye(size(temp_angle))))=inf;
                theta_min=min(temp_angle');
                theta_min=theta_min';
                
                arch = ArchiveUpdate(Population,Problem.N);
                Population = arch;
                gen  = 1;
                flag1 = 0;
                VAR0=0;
                    
                while Algorithm.NotTerminated(Population)
                     
                    Z = min(Population3.objs, [], 1); 
                    if flag1 == 0
                        std_obj(gen,:) = std(Population3.objs,[],1);
                        if gen>100
                            if  sum(std(std_obj(gen-100:gen,:),[],1)<0.1) == Problem.M
                                flag1 = 1;
                                Fes = Problem.FE;
                                cons = Population1.cons;
                                cons(cons<0) = 0;
                                cons =sum(cons,2);
                                index =find(cons>0);
                                if isempty(index)
                                    VAR0 = 0;
                                else
                                    VAR0 =  mean(cons(index));
                                end
                            end
                        end
                    end
               
                    for i = 1:Problem.N
                        ii = i;
                        if i > WVs
                            ii = randi([1,WVs]);
                        end
                        
                        [~,index1] = max(Q,[],2);
                        for j = 1:2
                            if isempty(find(Q(j,:)~=0, 1))
                                col1 = randperm(size(WVs,1));
                                index1(j) = col1(1);
                            end
                        end
    
                        if rand < greedy
                            PP{i} = randperm(WVs);
                            P = PP{i};
                            E(i,:)=P(1:2);
                            
                        else
                            PP{i} = B(ii,randperm(size(B,2)));
                            P = PP{i};
                            E(i,:)=P(1:2);
                        end
                        NN{i} = B(test_Q2(i),randperm(size(B,2)));
                        N = NN{i};
                    end
                    
    
                    Offspring1 = OperatorGA(Problem,Population1(randi(Problem.N,1,Problem.N)));
                    Offspring2 = OperatorGAhalf(Problem,[Population2(N(:,1)),Population2(N(:,2))]);
                    if flag1 == 0
    
                        Offspring3 = OperatorGAhalf(Problem,[Population3(E(:,1)),Population3(E(:,2))]); 
                        
                    else
                        Offspring3 = OperatorDE(Problem,Population3,Population3(E(:,1)),Population3(E(:,2)));
                    end      
                   Offspring = Offspring3;
    
                   
    
                   for i = 1 : Problem.N
                       
                       beta = 0.5 * (1 - cos((1 - gen / Problem.maxFE) * pi));
                       theta = 0.5*theta_min .* beta;
                       P = PP{i};
                       Z = min(Z,Offspring(i).obj);
                       PopObj1=Population3(P).objs-repmat(Z,length(P),1);
                       Angle0 = acos(1 - pdist2(real(PopObj1),W(P,:),'cosine'));
                       Angle2 = diag(Angle0);
                       NewAngle2=acos(1-pdist2(real(Offspring(i).objs-Z),W(P,:),'cosine'));
                       NewAngle2=NewAngle2';
                       
                       g_old2 = max(abs(Population3(P).objs-repmat(Z,length(P),1)).*W(P,:),[],2);
                       g_new2 = max(repmat(abs(Offspring(i).obj-Z),length(P),1).*W(P,:),[],2);
                       
                        
    
                       [~,index1] = max(Q,[],2);
                        for j = 1:2
                            if isempty(find(Q(j,:)~=0, 1))
                                col1 = randperm(size(WVs,1));
                                index1(j) = col1(1);
                            end
                        end
    
                        if rand > epsilon
                            if flag2 ==1
                                test_Q2(i) = index1(2);
    
                            else
                                test_Q2(i) = index1(1);
    
                            end
                        else
                            c = randperm(size(WVs,1));
                            test_Q2(i) = c(1);
                        end
                      
    
                       if (flag1 == 0)
                           case1=NewAngle2<= theta(P) & Angle2<= theta(P) & g_old2>=g_new2;
                           case2=NewAngle2> theta(P) & Angle2> theta(P) & g_old2>=g_new2;
                           case3=NewAngle2<= theta(P) & Angle2> theta(P)& g_old2>=g_new2;
                           case4=NewAngle2> theta(P) & Angle2<= theta(P) & g_old2>=g_new2;
                           Population3(P(find(case1 | case2 | case3 | case4 ,nr))) = Offspring(i);
                           if (find(case1 | case2 | case3 | case4 ,nr))
                               if flag2 ==1
                                    Q(2,test_Q2(i)) = (1-alp)*Q(2,test_Q2(i)) + alp*(2+gam*M(2,test_Q2(i))) + (1-alp)*gam*(max(Q(1,:))-M(2,test_Q2(i)));
                                    M(2,test_Q2(i)) = max(Q(1,:));
                               else
                                    Q(1,test_Q2(i)) = (1-alp)*Q(1,test_Q2(i)) + alp*(2+gam*M(1,test_Q2(i))) + (1-alp)*gam*(max(Q(1,:))-M(1,test_Q2(i)));
                                    M(1,test_Q2(i)) = max(Q(2,:));
                               end
                               flag2 = 1;
                               
                           else
                                if flag2 ==1
                                    Q(2,test_Q2(i)) = (1-alp)*Q(2,test_Q2(i)) + alp*(-0.5+gam*M(2,test_Q2(i))) + (1-alp)*gam*(max(Q(1,:))-M(2,test_Q2(i)));
                                    M(2,test_Q2(i)) = max(Q(1,:));
                                else
                                    Q(1,test_Q2(i)) = (1-alp)*Q(1,test_Q2(i)) + alp*(-0.5+gam*M(1,test_Q2(i))) + (1-alp)*gam*(max(Q(1,:))-M(1,test_Q2(i)));
                                    M(1,test_Q2(i)) = max(Q(2,:));
                                end
                                flag2 = 0; 
                           end
                               
                           
                       else
                           CVO2 = sum(max(0,Offspring(i).con));
                           CVP2 = sum(max(0,Population3(P).cons),2);
                           case1=NewAngle2<= theta(P) & Angle2<= theta(P) & (CVP2>CVO2 | (CVP2==CVO2 & g_old2>=g_new2));
                           case2=NewAngle2> theta(P) & Angle2> theta(P) & (CVP2>CVO2 | (CVP2==CVO2 & g_old2>=g_new2));
                           case3=NewAngle2<= theta(P) & Angle2> theta(P);
                           indices = mod(find([case1 ; case2 ; case3],nr),length(P));
                           if indices > 0
                               Population3(P(indices)) = Offspring(i);
                           end
                           if (find(case1 | case2 | case3 ,nr))
                               if flag2 ==1
                                    Q(2,test_Q2(i)) = (1-alp)*Q(2,test_Q2(i)) + alp*(2+gam*M(2,test_Q2(i))) + (1-alp)*gam*(max(Q(1,:))-M(2,test_Q2(i)));
                                    M(2,test_Q2(i)) = max(Q(1,:));
                               else
                                   Q(1,test_Q2(i)) = (1-alp)*Q(1,test_Q2(i)) + alp*(2+gam*M(1,test_Q2(i))) + (1-alp)*gam*(max(Q(1,:))-M(1,test_Q2(i)));
                                    M(1,test_Q2(i)) = max(Q(2,:));
                               end
                               flag2 = 1;
    
                           else
                                if flag2 ==1
                                    Q(2,test_Q2(i)) = (1-alp)*Q(2,test_Q2(i)) + alp*(-0.5+gam*M(2,test_Q2(i))) + (1-alp)*gam*(max(Q(1,:))-M(2,test_Q2(i)));
                                    M(2,test_Q2(i)) = max(Q(1,:));
                                else
                                    Q(1,test_Q2(i)) = (1-alp)*Q(1,test_Q2(i)) + alp*(-0.5+gam*M(1,test_Q2(i))) + (1-alp)*gam*(max(Q(1,:))-M(1,test_Q2(i)));
                                    M(1,test_Q2(i)) = max(Q(2,:));
                                end
                                flag2 = 0;
                           end
                       end
                   end
    
                 
                   [Population2,~] = EnvironmentalSelection([Population2,Offspring2,Offspring3],Problem.N,true);
                   if flag1 == 0
                       [Population1,~] = EnvironmentalSelection([Population1,Offspring1,Offspring2,Offspring3],Problem.N,true);
                       
                   else
                       [Population1,~] = EnvironmentalSelection_VAR([Population1,Offspring1,Offspring2,Offspring3,arch],Problem.N,VAR0); 
                       cons = Offspring1.cons;
                       cons(cons<0) = 0;
                       cons = sum(cons,2);
                       index = find(cons>0);
                       if isempty(index)
                           VAR0 = 0;
                       else
                           tanhFE = 8 * ((Problem.FE-Fes) / (Problem.maxFE-Fes)) - 4;
                           K = 0.5-tanh(tanhFE)/2;
                           VAR0 = K*mean(cons(index));
                       end
                   end
    
                   Population = [Population1,Population2,Population3];
                   arch = ArchiveUpdate([arch,Population],Problem.N);
               
                   gen =gen +1;
               
                   if Problem.FE >=Problem.maxFE
                       Population = arch;
                   end
                   
                end
            end
        end
    end