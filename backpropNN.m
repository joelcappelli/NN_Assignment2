% backprop a per-period backpropagation training for a multilayer feedforward
%          neural network.
%   Network = backprop(Layers,N,M,SatisfactoryMSE,Input,Desired) returns 
%   Network, a two field structure of the form Network.structure = Layers 
%   and Network.weights where weights is a cell array specifying the final 
%   weight matrices computed by minimizing the mean squared error between 
%   the Desired output and the actual output of the network given a set of 
%   training samples: Input and the SatisfactoryMSE (satisfactory mean 
%   squared error)
%
%   Input:
%    Layers - a vector of integers specifying the number of nodes at each
%     layer, i.e for all i, Layers(i) = number of nodes at layer i, there
%     must be at least three layers and the input layer Layers(1) must
%     equal the dimension of each vector in Input, likewise, Layers(end) 
%     must be equal to the dimension of each vector in Desired
%     N - training rate for network learning (0.1 - 0.9)
%     M - momentum for the weight update rule [0.1 - 0.9)
%     SatisfactoryMSE - the mse at which to terminate computation
%     Input - the training samples, a P-by-N matrix, where each Input[p] is
%      a training vector
%     Desired - the desired outputs, a P-by-M matrix where each Desired[p]
%      is the desired output for the corresponding input Input[p]
%
%   This algorithm uses the hyperbolic tangent node function 
%   2/(1+e^(-2*net)) - 1, suitable for use with bipolar data
%   
%   NOTE: due to its generality this algorithm is not as efficient as a 
%   one designed for a specific problem if the number of desired layers is 
%   known ahead of time, it is better to a) 'unfold' the loops inside the 
%   loop presenting the data. That is, calculate the input and output of each 
%   layer explicitly one by one and subsequently the modified error and weight 
%   matrix modifications b) remove momentum and training rate as parameters
%   if they are known

function Network = backpropNN(L,n,m,smse,X,D,fnHandles,maxNumEpochs,plotFig,hiddenBias, initWeights)

% determine number of input samples, desired output and their dimensions
[P,N] = size(X);
[Pd,M] = size(D);

ACT_FN = 1;
ACT_DERIV = 2;

% make user that each input vector has a corresponding desired output
if P ~= Pd 
    error('backprop:invalidTrainingAndDesired', ...
          'The number of input vectors and desired ouput do not match');
end

% make sure that at least 3 layers have been specified and that the 
% the dimensions of the specified input layer and output layer are
% equivalent to the dimensions of the input vectors and desired output
if length(L) < 3 
    error('backprop:invalidNetworkStructure','The network must have at least 3 layers');
else
    if N ~= L(1) || M ~= L(end)
        e = sprintf('Dimensions of input (%d) does not match input layer (%d)',N,L(1));
        error('backprop:invalidLayerSize', e);
    elseif M ~= L(end)
        e = sprintf('Dimensions of output (%d) does not match output layer (%d)',M,L(end));
        error('backprop:invalidLayerSize', e);    
    end
end

% will use the number of layers often, so save the number here
nLayers = length(L); 

% randomize the weight matrices (uniform random values in [-.5 .5], there
% is a weight matrix between each layer of nodes. Each layer (exclusive the 
% output layer) has a bias node whose activation is always 1, that is, the 
% node function is C(net) = 1. Furthermore, there is a link from each node
% in layer i to the bias node in layer j (the last row of each matrix)
% because this is less computationally expensive then the alternative.
% NOTE: below that the wieghts of all links to bias nodes are defined as
% zero
w = cell(nLayers-1,1); % a weight matrix between each layer
 
if(nargin >= 11)
    w = initWeights;
else
    for i=1:nLayers-2        
        %see Yann LeCun papers on setting weights 
        a = -sqrt(6/(L(i+1)+L(i)));
        b = -a;
        w{i} = [a + (b-a).*rand(L(i+1),L(i)+1) ; zeros(1,L(i)+1)];
    end

    %     a = -sqrt(6/(L(i+1)+L(i)));
    %     b = -a;
    %     w{i} = [a + (b-a).*rand(L(i+1),L(i)+1) 
    w{end} = 1 - 2.*rand(L(end),L(end-1)+1);
end

% preallocate activation,net vectors and delta weight matrices for faster 
% computation
% activation vectors, all but output layer include bias activation
a = cell(nLayers,1);
 
for i=1:nLayers-1
    a{i} = ones(L(i),1);
end
a{end} = ones(L(end),1);

% net vectors, one for each node in that layer but there is 
% no net for input layer
net = cell(nLayers-1,1);
for i=1:nLayers-2;
    net{i} = ones(L(i+1),1);
end
net{end} = ones(L(end),1);

% delta weight matrices
dw = cell(nLayers-1,1);
for i=1:nLayers-1
    dw{i} = zeros(size(w{i}));
end

% initialize stopping conditions
mse = Inf;  % assuming the intial weight matrices are bad
presentations = 0; % we'll measure by epoch instead of presentation
epochs = 0;

patternErrors = zeros(M,P*maxNumEpochs);
cycleErrors = zeros(1,maxNumEpochs);
cycleErr = 0;

% Uncomment to plot mse at each epoch
if(plotFig)
    fig1=figure;  
    mseArray = NaN(maxNumEpochs,1);
    epochArray = NaN(maxNumEpochs,1);
    h = plot(NaN,NaN,'r*');
    axis auto;
    grid on;
    xlabel('Epochs');
    ylabel('Mean square error')
    hold on;
end
while mse > smse && presentations < P * maxNumEpochs
    sse = 0; % running total of squared error
    %randomPresIndices = ceil(P*rand(P,1));
    for p=1:P 
        % get the current input vector and desired output
%         a{1} = X(randomPresIndices(p),:)';
%         Dp = D(randomPresIndices(p),:)';     
        a{1} = X(p,:)';
        Dp = D(p,:)';          
        % compute the inputs and outputs to each layer
        for i=1:nLayers-1
            % compute inputs to this layer
            net{i} = w{i} * a{i}; 

            % compute outputs of this layer
            % for all layers but the output layer, the last node is the 
            % bias node and its activation is 1
%             if i < nLayers-1 % inner layers
%                 a{i+1} = [G*tanh(A.*net{i}(1:end-1)) ; 1];            
%             else             % output layers
%                 a{i+1} = G*tanh(A.*net{i});            
%             end
            if i < nLayers-1 % inner layers
                    %a{i+1} = [fnHandles{ACT_FN}((net{i}(:,1:end))) hiddenBias*ones(P,1)];   
                    a{i+1} = [fnHandles{ACT_FN}(net{i}); hiddenBias];           
                else             % output layers
                    a{i+1} = fnHandles{ACT_FN}((net{i}));            
             end
        end
        
        % accumlate the squared error
        patternError = (Dp-a{end});
        sse = sse + sum(patternError.^2);
        
        % calculate the modified error at each layer and update the weight 
        % matrices accordingly. first calculate delta, the modified error
        % for the output nodes (S'(Output[net])*(Dp-Output[Activation])
        % then for each weight matrix, add n * delta * activation and
        % propagate delta to the previous layer
        %delta = (Dp-a{end}) .* (A/G).*(G-a{end}).*(G+a{end}); 
        delta = patternError.* fnHandles{ACT_DERIV}(a{end});  
        for i=nLayers-1:-1:1
            dw{i} = n * delta * a{i}' + (m .* dw{i});
            if i > 1 % dont compute mod err for input layer 
                %delta = (A/G).*(G-a{i}).*(G+a{i}).* (delta'*w{i})';
                delta = fnHandles{ACT_DERIV}(a{i}(1:end-1,:)).* (delta'*w{i}(:,1:end-1))';
                %delta = (A/G).*(G-a{i}).*(G+a{i}).* (delta'*w{i})';
            end
        end
            % update the prev_w, weight matrices, epoch count and mse
        for i=1:nLayers-1
            % we have the sum of the delta weights, divide through by the 
            % number of samples and add momentum * delta weight at (t-1)
            % finally, update the weight matrices
            w{i} = w{i} + dw{i};
        end     
        
        patternErrors(:,p + P*epochs) = 0.5*patternError.*patternError;
        cycleErr = cycleErr + sum(patternErrors(:,p + P*epochs),1);
    end
    presentations = presentations + P;
    mse = sse/(P*M); % mse = 1/P * 1/M * summed squared error
    epochs = epochs+1;
    
    cycleErrors(epochs) = cycleErr;
    cycleErr = 0;
    
    if(plotFig)
        mseArray(epochs)=mse;
        epochArray(epochs)=epochs;
        x=get(h,'XData');
        y=get(h,'YData');
        set(h,'XData',[x epochs],...
              'YData',[y mse]);
        drawnow
    end 
end

% return the trained network
Network.structure = L;
Network.weights = w;
Network.mse = mse;
Network.patternErrors = patternErrors;
Network.cycleErrors = cycleErrors;
Network.presentations = presentations;