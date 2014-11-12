data_set_1=[vector_1, vector_2];

data_set_1_sorted=sortrows(data_set_1);
data_set_1_sorted=flipdim(data_set_1_sorted,1);  

vector_1=data_set_1_sorted(:,1);
vector_2=data_set_1_sorted(:,2);