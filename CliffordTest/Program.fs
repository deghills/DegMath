open System
open DegMath
open Clifford

let testMultivectorPrint =
    Multivector.print

let testGetBlade : byte -> Multivector -> unit =
    fun b m -> Multivector.getBlade b m |> Console.WriteLine

let testVectors(p, q, n) =
    let cl = Cl(p, q, n)
    let vecs =
        cl.Vectors
        |> Seq.map (
            function KeyValue(potency, vecs) ->
                        vecs 
                        |> Seq.map (fun v -> potency, v)) 
                        |> Seq.collect id
    for (potency, v) in vecs do printfn $"{(potency, v)}"

let testBivectors(p, q, n) =
    for bv in Cl(p, q, n).Bivectors do printfn "%A" bv

let testDual : (int*int*int) -> Multivector -> unit =
    fun (p, q, n) ->
        Cl(p, q, n).Dual >> Multivector.print

let testDualInv : (int*int*int) -> Multivector -> unit =
    fun (p, q, n) ->
        Cl(p, q, n).DualInv >> Multivector.print

let testReverse : Multivector -> unit =
    Cl(3, 3, 2).Reverse >> Multivector.print

let testGradeProject : int -> Multivector -> unit =
    fun i m -> Cl(3, 3, 2).GradeProject i m |> Multivector.print

let testGrade : Multivector -> unit =
    Cl(3, 3, 2).Grade
    >> fun grades -> 
        for g in grades do Console.WriteLine g

let testNeg : Multivector -> unit =
    Cl(3, 3, 2).Neg >> Multivector.print

let testAdd : Multivector -> Multivector -> unit =
    fun a b -> Cl(3, 3, 2).Add a b |> Multivector.print

let testSub : Multivector -> Multivector -> unit =
    fun a b -> Cl(3, 3, 2).Sub a b |> Multivector.print

let testMul : (int*int*int) -> Multivector -> Multivector -> unit =
    fun (p, q, n) a b ->
        let cl = Cl(p, q, n)
        cl.Mul a b |> Multivector.print

let testDot : (int*int*int) -> Multivector -> Multivector -> unit =
    fun (p, q, n) a b ->
        let cl = Cl(p, q, n)
        cl.Dot a b |> Multivector.print

let testWdg : (int*int*int) -> Multivector -> Multivector -> unit =
    fun (p, q, n) a b ->
        let cl = Cl(p, q, n)
        cl.Wdg a b |> Multivector.print

let testReg : (int*int*int) -> Multivector -> Multivector -> unit =
    fun (p, q, n) a b ->
        let cl = Cl(p, q, n)
        cl.Reg a b |> Multivector.print

let testMagSqr : (int*int*int) -> Multivector -> unit =
    fun signature ->
        let cl = Cl signature
        cl.MagSqr >> Console.WriteLine

let testMag : (int*int*int) -> Multivector -> unit =
    fun signature ->
        let cl = Cl signature
        cl.Mag >> Console.WriteLine

let testNormalize : (int*int*int) -> Multivector -> unit =
    fun signature ->
        let cl = Cl signature
        cl.Normalize >> fst >> Multivector.print

let testVersorInv : Multivector -> unit =
    Cl(3, 3, 2).VersorInv >> Multivector.print

let testProject : (int*int*int) -> Multivector -> Multivector -> unit =
    fun signature a b ->
        Cl(signature).Project a b |> Multivector.print

let testSandwich : (int*int*int) -> Multivector -> Multivector -> unit =
    fun signature a b ->
        Cl(signature).Sandwich a b |> Multivector.print


let input1 = multivector [
    0b010uy, 1f
]

let input2 = multivector [
    0b011uy, 1f
    0b101uy, 1f
    0b110uy, 1f
]


[<EntryPoint>]
let main _ =
    let cl = (2, 0, 1)
    do
        Multivector.print input1 
        Console.WriteLine()

        testSandwich cl input1 input2

    0