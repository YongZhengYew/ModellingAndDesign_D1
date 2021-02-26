#GNU OCTAVE 6.1.0
1; # LET OCTAVE RECOGNIZE THIS IS NOT A FUNCTION FILE

function resIndex = getIndexAt(x, y, targetMatrix)
  len = rows(targetMatrix);
  resIndex = (x-1)*len + y;
  return;
endfunction

function resList = deconstructCoords(targetIndex, targetMatrix)
  len = rows(targetMatrix);
  resList = [floor(targetIndex/len)+1];
  resList = [resList, targetIndex-(resList(1)-1)*len];
endfunction

# create gaussian envelope
function resMatrix = SPDI(targetIndex,
                               targetMatrix,
                               heightFunc,
                               hillHeight)
  die = zeros(rows(targetMatrix), columns(targetMatrix));
  die(targetIndex) = 1;
  
  bwdie = bwdist(die);
  for i = 1:rows(bwdie)*columns(bwdie)
    targetMatrix(i) = targetMatrix(i) + heightFunc(bwdie(i), hillHeight);
  endfor
  
  resMatrix = targetMatrix;
endfunction

#callback for gaussHill
function res = cb_GAUSSIAN(currRadius, hillHeight)
  res = (1/normpdf(0))*hillHeight*normpdf((currRadius)/hillHeight);
endfunction

# dot product function
function resDot = windDotProduct(targetIndex,
                                 destinationIndex,
                                 targetMatrix,
                                 targetWindU,
                                 targetWindV)
  t = deconstructCoords(targetIndex, targetMatrix);
  d = deconstructCoords(destinationIndex, targetMatrix);
  dirx = d(1)-t(1);
  diry = d(2)-t(2);
  ls1 = [dirx, diry];
  ls2 = [targetWindU, targetWindV];
  resDot = dot(ls1, ls2)/5;
endfunction

# dot product with the whole map
function resMatrix = WDI(targetMatrix,
                                 targetWindUMap,
                                 targetWindVMap,
                                 destIndexList)
  resMatrix = zeros(rows(targetMatrix), columns(targetMatrix));
  for i = 1:rows(targetMatrix)*columns(targetMatrix)
    currSum = 0;
    count = 0;
    currWindU = targetWindUMap(i);
    currWindV = targetWindVMap(i);
    for destIndex = 1:length(destIndexList)
      currSum = currSum - windDotProduct(i,
                                         destIndexList(destIndex),
                                         targetMatrix,
                                         currWindU,
                                         currWindV);
      count++;
    endfor
    resMatrix(i) = targetMatrix(i) + currSum/count;
  endfor
endfunction

# remove water squares
function resDie = remvWater(targetMatrix, templateList)
  water = targetMatrix;
  waterCPtr = 1;
  for i = 1:length(templateList)
    currTempVal = templateList(i);
    if currTempVal != Inf
      water(getIndexAt(currTempVal, waterCPtr, water)) = Inf;
    else
      waterCPtr++;
    endif
  endfor
  resDie = water;
endfunction

# bounding function for desirability index tolerance
function resMatrix = boundSafeVals(targetMatrix, maxVal)
  for i = 1:rows(targetMatrix)*columns(targetMatrix)
    if targetMatrix(i) > maxVal
      targetMatrix(i) = Inf;
    endif
  endfor
  resMatrix = targetMatrix;
endfunction



# Global rows and column values
G_ROWS = 32;
G_COLUMNS = 32;

# Constructing Northeast to Southwest wind vector field
windvecu = zeros(G_ROWS, G_COLUMNS);
for i = 1:rows(windvecu)*columns(windvecu)
  windvecu(i) = -4.277778*sin(pi/4);
endfor
windvecv = zeros(G_ROWS, G_COLUMNS);
for i = 1:rows(windvecv)*columns(windvecv)
  windvecv(i) = 4.277778*cos(pi/4);
endfor

# Manually creating fluctuations
for i = 1:12
  for j = 1:6
    windvecu(j,i)+=1;
  endfor
endfor

for i = 1:12
  for j = 29:32
    windvecu(j,i)-=1;
  endfor
endfor

for i = 2:15
  for j = 3:5
    windvecu(j,i)-=1;
    windvecv(j,i)-=2;
  endfor
endfor

for i = 9:14
  for j = 18:23
    windvecu(j,i)-=1;
    windvecv(j,i)+=2;
  endfor
endfor

for i = 5:12
  for j = 23:29
    windvecu(j,i)-=1;
    windvecv(j,i)-=2;
  endfor
endfor

for i = 6:19
  for j = 1:9
    windvecu(j,i)+=2;
    windvecv(j,i)-=3;
  endfor
endfor

# List of water square locations, Inf is an instruction to remvWater() to proceed to the next row
waterLocList = [
  1:11, Inf, 1:9, Inf, 1:8, Inf, 1:8, Inf, 1:7, Inf, 1:6, Inf, 1:6, 9, Inf, 1:6, 8, 19:20, Inf,1:8, 14:15, 19, Inf,1:8, 14:18, Inf,1:7, 14:17, 22, Inf,  1:7, 15:18, 21, Inf,  1:6, 16:20, 29, Inf,  1:6, 16:19, 29, Inf,  1:6, 17:21, 29:30, Inf,  1:5, 17:22, 30:31, Inf,  1:5, 18:23, 31:32, Inf,  1:5, 18:24, 30:32, Inf,  1:5, 19:26, 29:32, Inf,  1:6, 18:32, Inf,  1:6, 19:32, Inf,  1:6, 19:30, Inf,  1:5, 19:29, Inf,  1:5, 19:27, Inf,  1:2, 16:26, Inf,  1:2, 16:22, Inf,  1, 17:21, Inf,  1, 16:21, Inf,  1, 17:21, Inf,  17:21, Inf,  1, 17:22, Inf,1, 17:22
]




####################################################
#EXECUTION
####################################################

# Initialize blank grid
x = zeros(G_ROWS, G_COLUMNS);

# Values of h for SPDI assigned
C_TUAS_SOUTH_INCIN_PLANT = getIndexAt(6, 17, x);
x = SPDI(C_TUAS_SOUTH_INCIN_PLANT, x, @cb_GAUSSIAN, 20);

C_TUAS_SOUTH_DESAL_PLANT = getIndexAt(7, 15, x);
x = SPDI(C_TUAS_SOUTH_DESAL_PLANT, x, @cb_GAUSSIAN, 10);

C_KEPPEL_SEGHERS_WTE_PLANT = getIndexAt(6, 16, x);
x = SPDI(C_KEPPEL_SEGHERS_WTE_PLANT, x, @cb_GAUSSIAN, 15);

C_TUAS_INCIN_PLANT = getIndexAt(13, 4, x);
x = SPDI(C_TUAS_INCIN_PLANT, x, @cb_GAUSSIAN, 20);

dotted = WDI(x,
             windvecu,
             windvecv,
             [C_TUAS_SOUTH_INCIN_PLANT,
              C_TUAS_SOUTH_DESAL_PLANT,
              C_KEPPEL_SEGHERS_WTE_PLANT,
              C_TUAS_INCIN_PLANT]);

#DI
divgAdded = dotted - 5*divergence(windvecu, windvecv);
  
dryLand = remvWater(divgAdded, waterLocList);

bounded = boundSafeVals(fillRes, 30);
