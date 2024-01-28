import seldom
import os
import sys
sys.path.append("..")
from dimcheck.dimcheck_core import *

class TestDimCheckSi(seldom.TestCase):
    @classmethod
    def setUpClass(self):
        self.si=si

    def tearDown(self):
        os.system('rm output/*')
    
    def test_si_property(self):
        self.assertEqual(si.is_pretty,False)
        self.assertEqual(si.unit_system,"si")
        self.assertNotEqual(si.quant_def_file.find("dimcheck/si.json"),-1)
        self.assertNotEqual(si.serialized_file.find("dimcheck/si.pickle"),-1)


    def test_si_save(self):
        # self.assertFalse(si.clean())
        self.assertEqual(si.hash,"")
        self.assertEqual(si.save(),True)
        self.assertNotEqual(si.hash,"")
        self.assertEqual(si.clean(),True)
        self.assertEqual(si.hash,"")

        expr_file = "."+os.sep+"output/Expressions.csv"
        self.assertTrue(si.save_all_expr(expr_file))
        self.assertTrue(os.path.exists(expr_file))

        quant_file = "."+os.sep+"output/Quantities.csv"
        self.assertTrue(si.save_all_quant(quant_file))
        self.assertTrue(os.path.exists(quant_file))

    def test_display(self):
        self.assertIsNone(si.all_expr())
        self.assertIsNone(si.all_quant())
        
    def test_si_is_pretty(self):
        is_pretty=si.is_pretty
        is_quant=si.is_quant
        try:
            si.is_pretty=True
            si.is_quant=False
            self.assertEqual(si.unit("[l]"), "m")
            self.assertEqual(si.unit("[F]"), "kg*m/s²")
            self.assertEqual(si.unit("[G]"), "m³/(kg*s²)")
            self.assertEqual(si.unit("[eps_0]"), "A²*s⁴/(kg*m³)")
            self.assertEqual(si.unit("[q]"), "A*s")
            self.assertEqual(si.unit("[rho]"), "A*s/m³")
            self.assertEqual(si.unit("[l]"), si.unit("[l]"))
            self.assertEqual(si.unit("[l]"), "m")
            self.assertEqual(si.unit("[F]"), "kg*m/s²")
            self.assertEqual(si.unit("[G]"), "m³/(kg*s²)")
            self.assertEqual(si.unit("[eps_0]"), "A²*s⁴/(kg*m³)")
            self.assertEqual(si.unit("[q]"), "A*s")
            self.assertEqual(si.unit("[rho]"), "A*s/m³")
            self.assertEqual(si.formula("[t]",["[l]","[g]"]),"[t] = [l]⁰·⁵*[g]⁻⁰·⁵")
            self.assertEqual(si.omit_quant(lhs="[E]",rhs="[k]",omit_quant=["sqrt([hbar v]**2)"]),"[E] = [k]*sqrt([hbar v]²)")
            self.assertEqual(si.quant("[V/I]"),['[R]', '[R_K]'])
            self.assertEqual(si.unit("sqrt[m^2 l^3]"),'kg*m³ʴ²')
            self.assertEqual(si.unit("[m]**(1/3)"),'kg¹ʴ³')
        except Exception as e:
            raise e
        finally:
            si.is_pretty=is_pretty
            si.is_quant=is_quant

    def test_si_unit(self):
        is_quant=si.is_quant
        try:
            si.is_quant=False
            self.assertEqual(si.unit("[l]"), "m")
            self.assertEqual(si.unit("[F]"), "kg*m/s**2")
            self.assertEqual(si.unit("[G]"), "m**3/(kg*s**2)")
            self.assertEqual(si.unit("[eps_0]"), "A**2*s**4/(kg*m**3)")
            self.assertEqual(si.unit("[q]"), "A*s")
            self.assertEqual(si.unit("[rho]"), "A*s/m**3")

            self.assertEqual(si.unit("[m][m]"),"kg**2")
            self.assertEqual(si.unit("[l/l/m][m]"),"1")
            self.assertEqual(si.unit("[m m / l g]"),"kg**2/s**2")
            self.assertEqual(si.unit("kg  m / s"),"kg*m/s")
            self.assertEqual(si.unit("kg**2 m / s"),"kg**2*m/s")
            self.assertEqual(si.unit("kg**2 m / s s"),"kg**2*m")
            self.assertEqual(si.unit("kg^(1/2) m / s"),"sqrt(kg)*m/s")
            self.assertEqual(si.unit("[hbar (v) ** (-3)* e ** (-2.0)/h]"),"s**1.0/(A**2.0*m**3)")
        except Exception as e:
            raise e
        finally:
            si.is_quant=is_quant
        

        

        

    def test_si_dimension(self):
        is_quant=si.is_quant
        try:
            si.is_quant=False
            self.assertEqual(si.dimension("[l]"), "L")
            self.assertEqual(si.dimension("[F]"), "M*L/T**2")
            self.assertEqual(si.dimension("[G]"), "L**3/(M*T**2)")
            self.assertEqual(si.dimension("[eps_0]"), "I**2*T**4/(M*L**3)")
            self.assertEqual(si.dimension("[q]"), "I*T")
            self.assertEqual(si.dimension("[rho]"), "I*T/L**3")

            self.assertEqual(si.dim("[l]"), "L")
            self.assertEqual(si.dim("[F]"), "M*L/T**2")
            self.assertEqual(si.dim("[G]"), "L**3/(M*T**2)")
            self.assertEqual(si.dim("[eps_0]"), "I**2*T**4/(M*L**3)")
            self.assertEqual(si.dim("[q]"), "I*T")
            self.assertEqual(si.dim("[rho]"), "I*T/L**3")
        except Exception as e:
            raise e
        finally:
            si.is_quant=is_quant
        

    def test_si_is_dc(self):
        is_quant=si.is_quant
        try:
            si.is_quant=False
            self.assertEqual(si.is_dc("[l]","m"), True)
            self.assertEqual(si.is_dc("[l]","kg"), False)
            self.assertEqual(si.is_dc("[l/l/m]","kg"),False)
            self.assertEqual(si.is_dc("[l/l*m]","kg"),True)
            self.assertEqual(si.is_dc("[F]","kg*m/s**2"), True)
            self.assertEqual(si.is_dc("[F]","[m*a]"), True)
            self.assertEqual(si.is_dc("[F]","[E/l]"), True)
            self.assertEqual(si.is_dc("[F]","[E]"), False)
            self.assertEqual(si.is_dc("[F]","[m]*[a]"), True)
            self.assertEqual(si.is_dc("[F]","[m]*[v/t]"), True)
            self.assertEqual(si.is_dc("[F]","[q*EE]"), True)
            self.assertEqual(si.is_dc("[F]","[q*E]"), False)
            self.assertEqual(si.is_dc("[q*EE]","[m*a]"), True)
            self.assertEqual(si.is_dc("[hbar v k]","[E]"),True)
            self.assertEqual(si.is_dc("[(hbar e B v**2)**(1/2)]","[E]"),True)
            self.assertEqual(si.is_dc("[hbar e B v**2]**0.5","[E]"),True)
            self.assertEqual(si.is_dc("[sqrt(hbar e B v**2)]","[E]"),True)
            self.assertEqual(si.is_dc("sqrt([v])**2","[v]"),True)
            self.assertEqual(si.is_dc("sqrt([v p E])","[E]"),True)
            self.assertEqual(si.is_dc("cbrt([v])**3","[v**(1/3)]**3"),True)
            self.assertEqual(si.is_dc("[mu_mob]","[l^2 /(V t)]",True),True)
            self.assertEqual(si.is_dc("[sigma]","[1/rho_R]",True),True)
            self.assertEqual(si.is_dc("[sigma_2D]","[e^2/hbar]",True),True)
            
            self.assertRaises(UnknownQuantityError,si.is_dc,"[q*EE]","[ma]")
        except Exception as e:
            raise e
        finally:
            si.is_quant=is_quant
        
        

    def test_si_quant(self):
        is_quant=si.is_quant
        try:
            si.is_quant=False
            self.assertEqual(si.quant("[l/t]"), ['[v]', '[c]'])
            self.assertEqual(si.quant("[r/t]"), ['[v]', '[c]'])
            self.assertEqual(si.quant("[r]"), "[l]")
            self.assertEqual(si.quant("[V/I]"), ['[R]', '[R_K]'])
            self.assertEqual(si.quant("[F*l]"), "[E]")
            self.assertEqual(si.quant("[G*m**2/r**2]"), "[F]")
            self.assertEqual(si.quant("[l/t**2]"), ['[a]', '[g]'])
            self.assertEqual(si.quant("[k_B*T]"), "[E]")
            self.assertEqual(si.quant("[E/q]"), "[V]")
            self.assertEqual(si.quant("[e]*[hbar]/[m]"), "[mu_B]")
            self.assertEqual(si.quant("kg m/ s"), "[p]")
        except Exception as e:
            raise e
        finally:
            si.is_quant=is_quant

    def test_si_omit_quant(self):
        is_quant=si.is_quant
        try:
            si.is_quant=False
            self.assertEqual(si.omit_quant("[E]","[v*k]","[hbar]"), "[E] = [v*k]*[hbar]")
            self.assertEqual(si.omit_quant("[E]**2","[v*k]**2","[hbar]"), "[E]**2 = [v*k]**2*[hbar]**2")
            self.assertEqual(si.omit_quant("[E]","[k]",["[hbar]","[c]"]), "[E] = [k]*[hbar]*[c]")
        except Exception as e:
            raise e
        finally:
            si.is_quant=is_quant

    def test_formula(self):
        is_quant=si.is_quant
        try:
            si.is_quant=False
            self.assertEqual(si.formula("[E]",["[v]","[k]","[hbar]"]), "[E] = [v]*[k]*[hbar]")
            self.assertEqual(si.formula("[E]",["[m]","[l]","[v]"]), "[E] = [v]**2*[m]")
            self.assertEqual(si.formula("[r]",["[G*m]","[t]"]), "[r] = [G*m]**(0.333333333)*[t]**(0.666666667)")
            self.assertEqual(si.formula("[t]",["[l]","[g]"]), "[t] = [l]**(0.5)*[g]**(-0.5)")
            self.assertEqual(si.formula("[G]",["[m]","[r]","[F]"]), "[G] = [m]**-2*[r]**2*[F]")
            self.assertEqual(si.formula("[t]",["[m/t**2]", "[V_l]", "[rho_m]","[g]"]),"[t] = [g]**(-0.5)*[V_l]**(0.166666667)")
            self.assertEqual(si.formula("[pressure]",["[m]","[n]","[T]","[m l^2 t^-2 T^-1]"]), "[pressure] = [m l^2 t^-2 T^-1]*[T]*[n]")
            self.assertEqual(si.formula("[pressure]",["[m]","[n]","[m l^2 t^-2]"]), "[pressure] = [m l^2 t^-2]*[n]")
            self.assertEqual(si.formula("[t]",["[I]","[L]","[C]"]), "[t] = [C]**(0.5)*[L]**(0.5)")
            self.assertEqual(si.formula("[t]",["[m]","[l]","[I]","[T]","[mole]","[F]","[a]","[Q]"]), "[t] = [Q]*[I]**-1")
            
        except Exception as e:
            raise e
        finally:
            si.is_quant=is_quant


    def test_si_is_quant(self):
        is_quant=si.is_quant
        is_pretty=si.is_pretty

        try:

            si.is_quant=True
            self.assertEqual(si.unit("l"), "m")

            self.assertEqual(si.unit("l"), "m")
            self.assertEqual(si.unit("F"), "kg*m/s**2")
            self.assertEqual(si.unit("G"), "m**3/(kg*s**2)")
            self.assertEqual(si.unit("eps_0"), "A**2*s**4/(kg*m**3)")
            self.assertEqual(si.unit("q"), "A*s")
            self.assertEqual(si.unit("rho"), "A*s/m**3")

            self.assertEqual(si.unit("m m"),"kg**2")
            self.assertEqual(si.unit("l/l/m m"),"1")
            self.assertEqual(si.unit("m m / l g"),"kg**2/s**2")
            self.assertEqual(si.unit("m  l / t"),"kg*m/s")
            self.assertEqual(si.unit("m**2 l / t"),"kg**2*m/s")
            self.assertEqual(si.unit("m**2 l / t t"),"kg**2*m")
            self.assertEqual(si.unit("m^(1/2) l / t"),"sqrt(kg)*m/s")
            self.assertEqual(si.unit("hbar (v) ** (-3)* e ** (-2.0)/h"),"s**1.0/(A**2.0*m**3)")

            self.assertEqual(si.is_dc("l","m"), False)
            self.assertEqual(si.is_dc("l/l/m","m"),False)
            self.assertEqual(si.is_dc("l/l*m","m"),True)
            self.assertEqual(si.is_dc("F","m*l/t**2"), True)
            self.assertEqual(si.is_dc("F","m*a"), True)
            self.assertEqual(si.is_dc("F","E/l"), True)
            self.assertEqual(si.is_dc("F","E"), False)
            self.assertEqual(si.is_dc("F","m*a"), True)
            self.assertEqual(si.is_dc("F","m*v/t"), True)
            self.assertEqual(si.is_dc("F","q*EE"), True)
            self.assertEqual(si.is_dc("F","q*E"), False)
            self.assertEqual(si.is_dc("q*EE","m*a"), True)
            self.assertEqual(si.is_dc("hbar v k","E"),True)
            self.assertEqual(si.is_dc("(hbar e B v**2)**(1/2)","E"),True)
            self.assertEqual(si.is_dc("(hbar e B (v**2))**0.5","E"),True)
            self.assertEqual(si.is_dc("sqrt(hbar e B v**2)","E"),True)
            self.assertEqual(si.is_dc("sqrt(v)**2","v"),True)
            self.assertEqual(si.is_dc("sqrt(v p E)","E"),True)
            self.assertEqual(si.is_dc("cbrt(v)**3","(v**(1/3))**3"),True)
            self.assertEqual(si.is_dc("mu_mob","l^2 /(V t)",True),True)
            self.assertEqual(si.is_dc("sigma","1/rho_R",True),True)
            self.assertEqual(si.is_dc("sigma_2D","e^2/hbar",True),True)
            
            self.assertRaises(UnknownQuantityError,si.is_dc,"q*EE","ma")

            self.assertEqual(si.quant("l/t"), ['v', 'c'])
            self.assertEqual(si.quant("r/t"), ['v', 'c'])
            self.assertEqual(si.quant("r"), "l")
            self.assertEqual(si.quant("V/I"), ['R', 'R_K'])
            self.assertEqual(si.quant("F*l"), "E")
            self.assertEqual(si.quant("G*m**2/r**2"), "F")
            self.assertEqual(si.quant("l/t**2"), ['a', 'g'])
            self.assertEqual(si.quant("k_B*T"), "E")
            self.assertEqual(si.quant("E/q"), "V")
            self.assertEqual(si.quant("e*hbar/m"), "mu_B")
            self.assertEqual(si.quant("m l/ t"), "p")

            self.assertEqual(si.omit_quant("E","v*k","hbar"), "E = v*k*hbar")
            self.assertEqual(si.omit_quant("E**2","(v*k)**2","hbar"), "E**2 = (v*k)**2*hbar**2")
            self.assertEqual(si.omit_quant("E","k",["hbar","c"]), "E = k*hbar*c")

            self.assertEqual(si.formula("E",["v","k","hbar"]), "E = v*k*hbar")
            self.assertEqual(si.formula("E",["m","l","v"]), "E = v**2*m")
            self.assertEqual(si.formula("r",["G m","t"]), "r = (G m)**(0.333333333)*t**(0.666666667)")
            self.assertEqual(si.formula("t",["l","g"]), "t = l**(0.5)*g**(-0.5)")
            self.assertEqual(si.formula("G",["m","r","F"]), "G = m**-2*r**2*F")
            self.assertEqual(si.formula("t",["m/t**2", "V_l", "rho_m","g"]),"t = g**(-0.5)*V_l**(0.166666667)")
            self.assertEqual(si.formula("pressure",["m","n","T","m l^2 t^-2 T^-1"]), "pressure = m l^2 t^-2 T^-1*T*n")
            self.assertEqual(si.formula("pressure",["m","n","m l^2 t^-2"]), "pressure = m l^2 t^-2*n")
            self.assertEqual(si.formula("t",["I","L","C"]), "t = C**(0.5)*L**(0.5)")
            self.assertEqual(si.formula("t",["m","l","I","T","mole","F","a","Q"]), "t = Q*I**-1")

            si.is_pretty=True
            self.assertEqual(si.unit("l"), "m")
            self.assertEqual(si.unit("F"), "kg*m/s²")
            self.assertEqual(si.unit("G"), "m³/(kg*s²)")
            self.assertEqual(si.unit("eps_0"), "A²*s⁴/(kg*m³)")
            self.assertEqual(si.unit("q"), "A*s")
            self.assertEqual(si.unit("rho"), "A*s/m³")
            self.assertEqual(si.unit("l"), "m")
            self.assertEqual(si.unit("F"), "kg*m/s²")
            self.assertEqual(si.unit("G"), "m³/(kg*s²)")
            self.assertEqual(si.unit("eps_0"), "A²*s⁴/(kg*m³)")
            self.assertEqual(si.unit("q"), "A*s")
            self.assertEqual(si.unit("rho"), "A*s/m³")
            self.assertEqual(si.formula("t",["l","g"]),"t = l⁰·⁵*g⁻⁰·⁵")
            self.assertEqual(si.omit_quant(lhs="E",rhs="k",omit_quant=["sqrt((hbar v)**2)"]),"E = k*sqrt((hbar v)²)")
            self.assertEqual(si.quant("V/I"),['R', 'R_K'])
            self.assertEqual(si.unit("sqrt(m^2 l^3)"),'kg*m³ʴ²')
            self.assertEqual(si.unit("m**(1/3)"),'kg¹ʴ³')
        except Exception as e:
            raise e
        finally:
            si.is_quant=is_quant
            si.is_pretty=is_pretty
        


    def test_si_error(self):
        is_quant=si.is_quant
        try:
            si.is_quant=False
            self.assertRaises(SquareBracketMismatchError,si.dim,"l]")
            self.assertRaises(SquareBracketMismatchError,si.dim,"m*a]/kg*[q]/[a")
            self.assertRaises(SquareBracketMismatchError,si.dim,"[l")
            self.assertRaises(MissingSquareBracketError,si.is_dc,"[q]","Q*[a]")
            self.assertRaises(MissingSquareBracketError,si.is_dc,"[q]","[Q]*Q")
            self.assertRaises(MissingQuantError,si.is_dc,"[q]/[]","[m]")
            self.assertRaises(MissingQuantError,si.is_dc,"[]","[m]")
            self.assertRaises(InvalidSymFormatError,si.is_dc,"[345l]","[m][a]")
            self.assertRaises(InvalidSyntaxError,si.is_dc,"[l$a]","[m][a]")
            self.assertRaises(UnknownQuantityError,si.is_dc,"[k]","[m]*[abc]")
            self.assertRaises(MissingParenthesesError,si.dim,"[v**1/3*m**2]")
            self.assertRaises(DerivingQuantityError,si.is_dc,"sqrt([v])","[v**0.5]")
            self.assertRaises(DerivingQuantityError,si.is_dc,"[v**0.5]","[v**(1/2)]")
            self.assertRaises(DerivingQuantityError,si.is_dc,"cbrt([v])","[v**(1/3)]")
    
        except Exception as e:
            raise e
        finally:
            si.is_quant=is_quant

        

    def test_custom(self):
        custom = Dimcheck(quant_def_file="./data/custom.json",is_save=False)

        is_pretty=custom.is_pretty
        is_quant=custom.is_quant

        try:
            custom.is_quant=False
            custom.is_pretty=False
            self.assertEqual(custom.is_pretty,False)
            self.assertEqual(custom.unit_system,"si")
            self.assertNotEqual(custom.quant_def_file.find("custom.json"),-1)
            self.assertNotEqual(custom.serialized_file.find("custom.pickle"),-1)


            self.assertFalse(custom.clean())
            self.assertEqual(custom.hash,"")
            self.assertEqual(custom.save(),True)
            self.assertNotEqual(custom.hash,"")
            self.assertEqual(custom.clean(),True)
            self.assertEqual(custom.hash,"")

            expr_file = "."+os.sep+"output/Expressions.csv"
            self.assertTrue(custom.save_all_expr(expr_file))
            self.assertTrue(os.path.exists(expr_file))

            quant_file = "."+os.sep+"output/Quantities.csv"
            self.assertTrue(custom.save_all_quant(quant_file))
            self.assertTrue(os.path.exists(quant_file))
            self.assertEqual(custom.unit("[l]"), "m")
            self.assertEqual(custom.unit("[F]"), "kg*m/s**2")
            self.assertEqual(custom.unit("[G]"), "m**3/(kg*s**2)")
            self.assertEqual(custom.unit("[eps_0]"), "A**2*s**4/(kg*m**3)")
            self.assertEqual(custom.unit("[q]"), "A*s")
            self.assertEqual(custom.unit("[rho]"), "A*s/m**3")
            

            


            self.assertEqual(custom.is_dc("[l]","m"), True)
            self.assertEqual(custom.is_dc("[l]","kg"), False)
            self.assertEqual(custom.is_dc("[F]","kg*m/s**2"), True)
            self.assertEqual(custom.is_dc("[F]","[m*a]"), True)
            self.assertEqual(custom.is_dc("[F]","[E/l]"), True)
            self.assertEqual(custom.is_dc("[F]","[E]"), False)
            self.assertEqual(custom.is_dc("[F]","[m]*[a]"), True)
            self.assertEqual(custom.is_dc("[F]","[m]*[v/t]"), True)
            self.assertEqual(custom.is_dc("[F]","[q*EE]"), True)
            self.assertEqual(custom.is_dc("[F]","[q*E]"), False)
            self.assertEqual(custom.is_dc("[q*EE]","[m*a]"), True)
            self.assertRaises(UnknownQuantityError,custom.is_dc,"[q*EE]","[ma]")

            self.assertEqual(custom.quant("[V/I]"),['[R]', '[R_K]'])
            self.assertEqual(custom.quant("[F*l]"), "[E]")
            self.assertEqual(custom.quant("[G*m**2/r**2]"), "[F]")
            self.assertEqual(custom.quant("[l/t**2]"),['[a]', '[g]'])

            self.assertEqual(custom.formula("[r]",["[kappa]","[t]"]), "[r] = [kappa]**(0.333333333)*[t]**(0.666666667)")


            custom.is_pretty=True
            self.assertEqual(custom.unit("[l]"), "m")
            self.assertNotEqual(custom.quant("[m a^2]"), "[F]")
            self.assertEqual(custom.quant("[m a]"), "[F]")
            self.assertEqual(custom.unit("[F]"), "kg*m/s²")
            self.assertEqual(custom.unit("[G]"), "m³/(kg*s²)")
            self.assertEqual(custom.unit("[eps_0]"), "A²*s⁴/(kg*m³)")
            self.assertEqual(custom.unit("[q]"), "A*s")
            self.assertEqual(custom.unit("[rho]"), "A*s/m³")
            custom.is_pretty=False
        
        except Exception as e:
            raise e
        finally:
            custom.is_pretty=is_pretty
            custom.is_quant=is_quant


    def test_custom_error(self):
        ## multi-defined symbol [P] in lhs
        self.assertRaises(ConflictSymbolDefinitionError,Dimcheck,"./data/custom_ConflictSymbolDefinitionError_I.json",is_save=False)

        ## multi-defined symbol [P] in alias
        self.assertRaises(ConflictSymbolDefinitionError,Dimcheck,"./data/custom_ConflictSymbolDefinitionError_II.json",is_save=False)

        ## [mkl] at the middle of the file cannot be derived
        self.assertRaises(DerivingQuantityError,Dimcheck,"./data/custom_DerivingQuantityError_I.json",is_save=False)

        ## [ml] at the end of the file cannot be derived
        self.assertRaises(DerivingQuantityError,Dimcheck,"./data/custom_DerivingQuantityError_II.json",is_save=False)




if __name__ == "__main__":
    seldom.main(open=False)

        