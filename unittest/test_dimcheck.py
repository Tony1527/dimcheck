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
        si.is_pretty=True
        self.assertEqual(si.dim("[l]"), "m")
        self.assertEqual(si.dim("[F]"), "kg*m/s²")
        self.assertEqual(si.dim("[G]"), "m³/(kg*s²)")
        self.assertEqual(si.dim("[eps_0]"), "A²*s⁴/(kg*m³)")
        self.assertEqual(si.dim("[q]"), "A*s")
        self.assertEqual(si.dim("[rho]"), "A*s/m³")
        self.assertEqual(si.dimension("[l]"), si.dim("[l]"))
        self.assertEqual(si.dimension("[l]"), "m")
        self.assertEqual(si.dimension("[F]"), "kg*m/s²")
        self.assertEqual(si.dimension("[G]"), "m³/(kg*s²)")
        self.assertEqual(si.dimension("[eps_0]"), "A²*s⁴/(kg*m³)")
        self.assertEqual(si.dimension("[q]"), "A*s")
        self.assertEqual(si.dimension("[rho]"), "A*s/m³")
        self.assertEqual(si.formula("[t]",["[l]","[g]"]),"[t] = ([l]⁰·⁵)*([g]⁻⁰·⁵)")
        self.assertEqual(si.omit_quant(lhs="[E]",rhs="[k]",omit_quant=["sqrt([hbar v]**2)"]),"[E] = [k]*sqrt([hbar v]²)")
        self.assertEqual(si.quant("[V/I]"),['[R]', '[R_K]'])
        self.assertEqual(si.dim("sqrt[m^2 l^3]"),'kg*m³ʴ²')
        self.assertEqual(si.dim("[m]**(1/3)"),'kg¹ʴ³')
        si.is_pretty=False

    def test_si_dim(self):
        self.assertEqual(si.dim("[l]"), "m")
        self.assertEqual(si.dim("[F]"), "kg*m/s**2")
        self.assertEqual(si.dim("[G]"), "m**3/(kg*s**2)")
        self.assertEqual(si.dim("[eps_0]"), "A**2*s**4/(kg*m**3)")
        self.assertEqual(si.dim("[q]"), "A*s")
        self.assertEqual(si.dim("[rho]"), "A*s/m**3")

        self.assertEqual(si.dim("[m][m]"),"kg**2")
        self.assertEqual(si.dim("[l/l/m][m]"),"1")
        self.assertEqual(si.dim("[m m / l g]"),"kg**2/s**2")
        self.assertEqual(si.dim("kg  m / s"),"kg*m/s")
        self.assertEqual(si.dim("kg**2 m / s"),"kg**2*m/s")
        self.assertEqual(si.dim("kg**2 m / s s"),"kg**2*m")
        self.assertEqual(si.dim("kg^(1/2) m / s"),"sqrt(kg)*m/s")
        self.assertEqual(si.dim("[hbar (v) ** (-3)* e ** (-2.0)/h]"),"s**1.0/(A**2.0*m**3)")
        

        

        

    def test_si_dimension(self):
        self.assertEqual(si.dimension("[l]"), si.dim("[l]"))
        self.assertEqual(si.dimension("[l]"), "m")
        self.assertEqual(si.dimension("[F]"), "kg*m/s**2")
        self.assertEqual(si.dimension("[G]"), "m**3/(kg*s**2)")
        self.assertEqual(si.dimension("[eps_0]"), "A**2*s**4/(kg*m**3)")
        self.assertEqual(si.dimension("[q]"), "A*s")
        self.assertEqual(si.dimension("[rho]"), "A*s/m**3")
        

    def test_si_is_dc(self):
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

    def test_si_quant(self):
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

    def test_si_omit_quant(self):
        self.assertEqual(si.omit_quant("[E]","[v*k]","[hbar]"), "[E] = [v*k]*[hbar]")
        self.assertEqual(si.omit_quant("[E]**2","[v*k]**2","[hbar]"), "[E]**2 = [v*k]**2*([hbar]**2)")
        self.assertEqual(si.omit_quant("[E]","[k]",["[hbar]","[c]"]), "[E] = [k]*[hbar]*[c]")

    def test_formula(self):
        self.assertEqual(si.formula("[E]",["[v]","[k]","[hbar]"]), "[E] = [v]*[k]*[hbar]")
        self.assertEqual(si.formula("[E]",["[m]","[l]","[v]"]), "[E] = [m]*([v]**2)")
        self.assertEqual(si.formula("[r]",["[G*m]","[t]"]), "[r] = ([G*m]**0.333333333)*([t]**0.666666667)")
        self.assertEqual(si.formula("[t]",["[l]","[g]"]), "[t] = ([l]**0.5)*([g]**-0.5)")
        self.assertEqual(si.formula("[G]",["[m]","[r]","[F]"]), "[G] = ([m]**-2)*([r]**2)*[F]")


    def test_si_is_quant(self):
        si.is_quant=True
        self.assertEqual(si.dim("l"), "m")

        self.assertEqual(si.dim("l"), "m")
        self.assertEqual(si.dim("F"), "kg*m/s**2")
        self.assertEqual(si.dim("G"), "m**3/(kg*s**2)")
        self.assertEqual(si.dim("eps_0"), "A**2*s**4/(kg*m**3)")
        self.assertEqual(si.dim("q"), "A*s")
        self.assertEqual(si.dim("rho"), "A*s/m**3")

        self.assertEqual(si.dim("m m"),"kg**2")
        self.assertEqual(si.dim("l/l/m m"),"1")
        self.assertEqual(si.dim("m m / l g"),"kg**2/s**2")
        self.assertEqual(si.dim("m  l / t"),"kg*m/s")
        self.assertEqual(si.dim("m**2 l / t"),"kg**2*m/s")
        self.assertEqual(si.dim("m**2 l / t t"),"kg**2*m")
        self.assertEqual(si.dim("m^(1/2) l / t"),"sqrt(kg)*m/s")
        self.assertEqual(si.dim("hbar (v) ** (-3)* e ** (-2.0)/h"),"s**1.0/(A**2.0*m**3)")

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
        self.assertEqual(si.omit_quant("E**2","(v*k)**2","hbar"), "E**2 = (v*k)**2*(hbar**2)")
        self.assertEqual(si.omit_quant("E","k",["hbar","c"]), "E = k*hbar*c")

        self.assertEqual(si.formula("E",["v","k","hbar"]), "E = v*k*hbar")
        self.assertEqual(si.formula("E",["m","l","v"]), "E = m*(v**2)")
        self.assertEqual(si.formula("r",["(G*m)","t"]), "r = ((G*m)**0.333333333)*(t**0.666666667)")
        self.assertEqual(si.formula("t",["l","g"]), "t = (l**0.5)*(g**-0.5)")
        self.assertEqual(si.formula("G",["m","r","F"]), "G = (m**-2)*(r**2)*F")


        si.is_pretty=True
        self.assertEqual(si.dim("l"), "m")
        self.assertEqual(si.dim("F"), "kg*m/s²")
        self.assertEqual(si.dim("G"), "m³/(kg*s²)")
        self.assertEqual(si.dim("eps_0"), "A²*s⁴/(kg*m³)")
        self.assertEqual(si.dim("q"), "A*s")
        self.assertEqual(si.dim("rho"), "A*s/m³")
        self.assertEqual(si.dimension("l"), si.dim("l"))
        self.assertEqual(si.dimension("l"), "m")
        self.assertEqual(si.dimension("F"), "kg*m/s²")
        self.assertEqual(si.dimension("G"), "m³/(kg*s²)")
        self.assertEqual(si.dimension("eps_0"), "A²*s⁴/(kg*m³)")
        self.assertEqual(si.dimension("q"), "A*s")
        self.assertEqual(si.dimension("rho"), "A*s/m³")
        self.assertEqual(si.formula("t",["l","g"]),"t = (l⁰·⁵)*(g⁻⁰·⁵)")
        self.assertEqual(si.omit_quant(lhs="E",rhs="k",omit_quant=["sqrt((hbar v)**2)"]),"E = k*sqrt((hbar v)²)")
        self.assertEqual(si.quant("V/I"),['R', 'R_K'])
        self.assertEqual(si.dim("sqrt(m^2 l^3)"),'kg*m³ʴ²')
        self.assertEqual(si.dim("m**(1/3)"),'kg¹ʴ³')
        si.is_pretty=False

        si.is_quant=False


    def test_si_error(self):
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

        

    def test_custom(self):
        custom = Dimcheck(quant_def_file="./data/custom.json",is_save=False)
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
        self.assertEqual(custom.dim("[l]"), "m")
        self.assertEqual(custom.dim("[F]"), "kg*m/s**2")
        self.assertEqual(custom.dim("[G]"), "m**3/(kg*s**2)")
        self.assertEqual(custom.dim("[eps_0]"), "A**2*s**4/(kg*m**3)")
        self.assertEqual(custom.dim("[q]"), "A*s")
        self.assertEqual(custom.dim("[rho]"), "A*s/m**3")
        

        


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

        self.assertEqual(custom.formula("[r]",["[kappa]","[t]"]), "[r] = ([kappa]**0.333333333)*([t]**0.666666667)")


        custom.is_pretty=True
        self.assertEqual(custom.dim("[l]"), "m")
        self.assertNotEqual(custom.quant("[m a^2]"), "[F]")
        self.assertEqual(custom.quant("[m a]"), "[F]")
        self.assertEqual(custom.dim("[F]"), "kg*m/s²")
        self.assertEqual(custom.dim("[G]"), "m³/(kg*s²)")
        self.assertEqual(custom.dim("[eps_0]"), "A²*s⁴/(kg*m³)")
        self.assertEqual(custom.dim("[q]"), "A*s")
        self.assertEqual(custom.dim("[rho]"), "A*s/m³")
        custom.is_pretty=False


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

        