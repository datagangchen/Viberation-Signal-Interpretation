function st = get_params(phi)
  
  st = struct;
  
  switch (phi.type)   
   case 'predicate'
    st = rmfield(phi.params,'fn');
    
   case 'not'
    st = get_params(phi.phi);   
   
   case 'or'
    st1 =  get_params(phi.phi1);    
    st2 =  get_params(phi.phi2);        
    st = ConcatStruct( st1,st2);
    
   case 'and'
    st1 =   get_params(phi.phi1);    
    st2 =   get_params(phi.phi2);        
    st = ConcatStruct( st1,st2);
    
   case '=>'
    st1 =   get_params(phi.phi1);    
    st2 =   get_params(phi.phi2);        
    st = ConcatStruct( st1,st2);
    
   case 'always'
    st = get_params(phi.phi);    
   
   case 'eventually'
    st = get_params(phi.phi);   
    
   case 'until'
    st1 = get_params(phi.phi1);   
    st2 = get_params(phi.phi2);   
    st = ConcatStruct(st1,st1);
    
  end

function st1 = ConcatStruct(st1,st2)
  listfield2 = fieldnames(st2);
  fieldvalues2 = struct2cell(st2);
  
  for i = 1:numel(listfield2)
    st1 = setfield(st1, listfield2{i}, fieldvalues2{i});
  end

    