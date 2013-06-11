class LPTutSolver(object):
    _initialtab = 0
    _m = 0
    _n = 0
    _initialb = 0
    _initialc = 0
    _tut = False
    _numVariables = 0
    _variables = 0

    def __init__(self, variables, ob, con, maximize = True, tut = False):
        LPtab = self._lptab(variables, ob, con, maximize, tut)
        self._initialb = LPtab["b"]
        self._initialtab = LPtab["lptab"]
        self._initialc = LPtab["c"]
        self._m = self._initialtab.nrows()
        self._n = self._initialtab.ncols()
        self._tut = tut
        self._numVariables = len(variables)
        self._variables = variables
        print "Initial:"
        print self._initialtab
    ######################################
    def _lptab(self, variables, ob, con, maximize = True, tut = False):
        conlist = []
        A = [[] for i in range(len(con))]
        b = [0 for i in range(len(con))]
        if maximize:
            c = [ob.coeff(i) for i in variables]
        else:
            c = [-ob.coeff(i) for i in variables]
        for i in range(len(con)):
            opera = con[i].operator()
            if opera == operator.le:
                A[i] = [con[i].left().coeff(j) for j in variables]
                b[i] = con[i].right()
            elif opera == operator.ge:
                A[i] = [-(con[i].left().coeff(j)) for j in variables]
                b[i] = -con[i].right()
            else:
                A[i] = [con[i].left().coeff(j) for j in variables]
                b[i] = con[i].right()
                conlist.append(i)
        if len(conlist) > 0:
            for i in conlist:
                A.append([])
                A[len(b)] = [-(con[i].left().coeff(j)) for j in variables]
                b.append(-con[i].right())
        A = matrix(A)
        lptab = block_matrix([[A.augment(identity_matrix(A.nrows())),matrix(b).transpose()],[matrix(c).augment(zero_matrix(SR,1,A.nrows())),0]])
        if A.nrows() == 2 and A.ncolumns == 2:
            _plotter(A,b)
        return {"lptab": lptab, "A": A, "b": b, "c": c}
    ######################################
    def _plotter(self,A,b):       #This plots the constraints and objective only for 2D graphs
        intcpt1=self.b[0]/self.A[0][1]
        intcpt2=self.b[1]/self.A[1][1]
        slope1=self.A[0][0]/self.A[0][1]
        slope2=self.A[1][0]/self.A[1][1]
        f1=intcpt1-slope1*x
        f2=intcpt2-slope2*x
        root1=find_root(f1==f2,-len(b)^10,len(b)^10)
        root2=find_root(f1==0,-len(b)^10,len(b)^10)
        root3=find_root(f2==0,-len(b)^10,len(b)^10)
        if slope1 > slope2:
            if root1 > 0:
                plot(f1,(x,root2,root1),rgbcolor=hue(1.0),fill=true) + plot(f2,(x,root1,root3),fill=true)
            elif root1==0:
                plot(f1,(x,root1,abs(max(b))),rgbcolor=hue(1.0),fill=true) + plot(f2, (x,root1,abs(max(b))),fill=true)
            else:
                plot(f1,(x,root2,root1),rgbcolor=hue(1.0),fill=true) + plot(f2, (x,root3,root1),fill=true)
        elif slope1 < slope2:
            if root1 > 0:
                plot(f1,(x,root1,abs(max(b))),rgbcolor=hue(1.0),fill=true) + plot(f2,(x,root3,root1), fill=true)
            elif root1==0:
                plot(f1,(x,root1,abs(max(b))),rgbcolor=hue(1.0),fill=true) + plot(f2, (x,root1,abs(max(b))),fill=true)
            else:
                plot(f1,(x,root2,root1),rgbcolor=hue(1.0),fill=true) + plot(f2, (x,root3,root1),fill=true)
    ######################################
    def _testers(self, b, c):
        ####Test to see if Tableaux is Primal Feasible
        testP = True
        for i in b:
            if i < 0: #if b is all positive we can proceed with Primal simplex
                testP = False
                break
        ####Test to see if tableaux is Dual Feasible
        testD = False
        if testP != True:
            testD = True
            for i in c:
                if i >= 0: #if last row is all negative we can proceed with Dual Simplex
                    testD = False
                    break

        if testP == False and testD == False:
            print "The linear program is infeasible. We are done."
        return (testP, testD)
    #######################################
    def _pivcol(self,tableaux):
        #choose pivot column based on Bland's rule of choosing left most positive value of c
        for i in range(0,self._n):
            if tableaux[self._m-1, i] > 0:
                pivotcolindex = i
                break
            else:
                i += 1
        return pivotcolindex
    #######################################
    def _pivrow(self, b, tableaux, pivotcolindex):
        #choose pivot row based on ratios
        pivotcol = tableaux.column(pivotcolindex).list()
        ratios = [] #calculate ratios
        a = "b/Zero"
        for i in range(len(b)): #only calculate ratios for non-objective rows
            if pivotcol[i] > 0: #makes sure we find ratio that is not dividing by 0 or dividing by a negative.
                ratio = b[i]/pivotcol[i]
                ratios.append(ratio)
            else:
                ratios.append(a)
        pivotrowindex = ratios.index(min(ratios)) #find index of the minimum ratio. This is our row
        if self._tut:
            print "index of pivotrow: " + str(pivotrowindex)
            print "index of pivotcol: " + str(pivotcolindex)
            print "ratios are: " + str(ratios)

        return pivotrowindex
    #######################################
    def _stepG(self,pivotrowindex,pivotcolindex, tableaux):
        #Now that we have the pivot row and pivot column, we can implement Gaussian elimination.
        #Find a, alpha, b
        pivotcol = tableaux.column(pivotcolindex).list()
        a = pivotcol[:pivotrowindex] #define values above pivot point, alpha
        alpha = pivotcol[pivotrowindex] #define alpha
        b = pivotcol[pivotrowindex+1:] #define values below pivot point, alpha
        if self._tut:
            print "a: " + str(a)
            print "alpha: " + str(alpha)
            print "b: " + str(b)
        #Now Build G
        #building top row of G
        aM = Matrix(SR,len(a), a) #turns a into a matrix for future augmenting
        if len(b) != 0: #created if statement to account for when b is the empty vector because changes the structure of G
            I_s = identity_matrix(SR, len(a)) #creates matrix of zeros
            alphaa = (-alpha**-1)*aM
            szeros = zero_matrix(SR, len(a), len(b))
            toprow = I_s.augment(alphaa.augment(szeros))
        else:
            I_s = identity_matrix(SR, len(a)) #creates matrix of zeros
            alphaa = (-alpha**-1)*aM
            toprow = I_s.augment(alphaa)
        #building middle row of G
        center = [alpha**-1]
        centerM = Matrix(SR, 1, 1, [center])
        middlezerosL = zero_matrix(SR, 1, len(a))
        middlezerosR = zero_matrix(SR, 1, len(b))
        middlerow = middlezerosL.augment(centerM.augment(middlezerosR))
        #building bottom row of G
        bM = Matrix(SR, b).transpose() #changes b from a list to a matrix to augment into bottomrow
        if len(a) != 0: #created if statement to account for when a is the empty vector because changes G's structure
            tzeros = zero_matrix(SR, len(b), len(a)) #build zero matrix
            alphab = (-alpha**-1)*bM #builds alpha matrix
            I_t = identity_matrix(SR, len(b)) #builds identity matrix (t X t)
            bottomrow = tzeros.augment(alphab.augment(I_t)) #augments to build bottom row of G
        else:
            alphab = (-alpha**-1)*bM #builds alpha matrix
            I_t = identity_matrix(SR, len(b)) #builds identity matrix (t X t)
            bottomrow = alphab.augment(I_t) #augments to build bottom row of G
        if len(a) == 0: #Changes the structure of G depending if a or b are the empty vector
            G = middlerow.stack(bottomrow)
        elif len(b)== 0:
            G = toprow.stack(middlerow)
        else:
            G = toprow.stack(middlerow.stack(bottomrow))
        if self._tut:
            print "Gaussian elimination matrix: "
            print G
        return G
    #######################################
    def _pivot(self, G, tableaux):
        #Now multiply Gaussian matrix by [a, alpha, b] to get the pivoted tableau
        return G * tableaux
    #######################################
    def _dualpivrow(self, tableaux):
        #choose pivot row based on largest negative of constraint column
        pivotrowindex = 0
        largestnegative = 0 #just put an arbitrarily high number to make the checks work. Assumes all constraints are less than 10^10
        for i in range(0,self._m):
            if tableaux[i, self._n-1] < 0:
                if tableaux[i, self._n-1] < largestnegative:
                       largestnegative = tableaux[i, self._n-1]
                       pivotrowindex = i
        return pivotrowindex #index of pivot row
    #######################################
    def _dualpivcol(self, tableaux, pivotrowindex):
        #Now that we have pivot row, find pivot column based on ratios.
        pivotrow = list(tableaux[pivotrowindex])
        objectiverow = list(tableaux[self._m-1])
        print pivotrow
        print objectiverow
        #Note, this ratio portion needs some help on eliminating ratios that are not allowed
        ratios = [] #calculate ratios
        for i in range(self._n - 1): #only calculate ratios for non-objective rows
            if pivotrow[i] != 0:
                ratio = objectiverow[i]/pivotrow[i]
                ratios.append(-1 * ratio)
            else:
                ratios.append("b/Zero")
        pivotcolindex = ratios.index(min(ratios)) #find index of the minimum ratio. This is our row
        if self._tut:
            print "index of pivotrow: " + str(pivotrowindex)
            print "index of pivotcol: " + str(pivotcolindex)
            print "ratios are: " + str(ratios)
        return pivotcolindex
    #######################################
    def solver(self):
        tableaux = self._initialtab
        b = self._initialb
        c = self._initialc
        testres = self._testers(b,c)
        testP = testres[0]
        testD = testres[1]
        counter = 1
        if testP:                     #preforming tests in the primal tableaux
            testdone = False
            while testdone == False:
                if self._tut:
                    print "Iteration #" + str(counter)
                pivotcol = self._pivcol(tableaux)
                pivotrow = self._pivrow(b, tableaux, pivotcol)
                G = self._stepG(pivotrow, pivotcol, tableaux)
                tableaux = self._pivot(G, tableaux)
                b = tableaux.column(self._n-1).list()[:self._m-1]
                testdone = True
                for i in tableaux.rows()[self._m-1]: #### testing to see if we are finished or not
                    if i > 0:
                        testdone = False
                        break
                if testdone == True:
                    break
                elif counter > 25:
                    break
                else:                                #if it is not done, print off the feasible soltn/objective
                    if self._tut:
                        line = "A feasible solution of the Primal is:\n"
                        line = line + self._printer(tableaux)
                        objectivevalue = -1*tableaux[self._m-1][self._n-1]
                        line = line + "The objective value of the feasible solution is:\n"
                        print line + str(objectivevalue)
                        print ""
                    counter += 1
        elif testD:            #preforming tests in the dual tableaux
            testdone = False
            while testdone == False:
                if self._tut:
                    print "Iteration #" + str(counter)
                pivotrow = self._dualpivrow(tableaux)
                pivotcol = self._dualpivcol(tableaux, pivotrow)
                G = self._stepG(pivotrow, pivotcol, tab)
                tableaux = self._pivot(G, tableaux)
                b = tableaux.column(self._n-1).list()[:self._m-1]
                testdone = True
                for i in tableaux.columns()[self._m-1]: #### testing to see if we are finished or not
                    if i < 0:
                        testdone = False
                        break
                if testdone == True:
                    break
                elif counter > 25:
                    break
                else:                    #prints the dual fesaible soltn/value
                    if self._tut:
                        line = "A feasible solution of the Dual is:\n"
                        line = line + self._printer(tableaux) + "\n"
                        line = line + "The objective value of the feasible solution is:\n"
                        print line + str(tableaux[(self._m-1),(self._n-1)])
                        print ""
                    counter += 1
        line = "The optimal solution to the Primal is:\n"        #once the algorithm is finished we print final soltn and max/min value
        line = line + self._printer(tableaux)
        line = line + "The maximum objective value of the Primal is:\n"
        m = -1
        if testD:
            m = 1
        objectivevalue = m*tableaux[self._m-1][self._n-1]
        print tableaux
        print line + str(objectivevalue)
    #######################################
    def _printer(self, tableaux):
        output = ""
        for i in range(0,self._numVariables):
            if tableaux[self._m-1][i] == 0: #identify if the columns of the variables have 0 in the objective row
                basisxcol = tableaux[:,i] #extract the column corresponding to x variables in the basis.
                for j in range(0, self._m-1): #loop over basis column to identify row index of the identity, 1
                    if basisxcol[j][0] == 1:
                        output = output + str(self._variables[i]) + " = " + str(tableaux[j][self._n-1]) + "\n" #find the corresponding value in the constraint column
            else:
                output = output + str(self._variables[i]) + " = 0\n"
        return output
    #######################################
    def toggleTut(self):
        self._tut = not self._tut
