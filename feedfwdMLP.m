% feed forward multi-layer perceptron

function output = feedfwdMLP(X,L,w,fnHandles,hiddenBias)

% determine number of input samples, desired output and their dimensions
[P,N] = size(X);
output = zeros(P,L(end));

ACT_FN = 1;
ACT_DERIV = 2;

% will use the number of layers often, so save the number here
nLayers = length(L); 

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

for p=1:P 
     a{1} = X(p,:)';
    % compute the inputs and outputs to each layer
    for i=1:nLayers-1
        % compute inputs to this layer
        net{i} = w{i} * a{i}; 

        if i < nLayers-1 % inner layers
                a{i+1} = [fnHandles{ACT_FN}(net{i}); hiddenBias];           
            else             % output layers
                a{i+1} = fnHandles{ACT_FN}((net{i}));            
         end
    end

    % accumlate the squared error
    output(p,:) = a{end};
end
end

