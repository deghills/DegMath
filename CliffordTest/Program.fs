open DegMath
open Clifford
open Multivector

let testMultivectorPrint(m:Multivector) =
    Multivector.print m

let testMultivectorOfBlade =
    Multivector.ofBlade >> Multivector.print

let testMultivectorOfBlades : Blade seq -> unit =
    Multivector.ofBlades >> Multivector.print

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


[<EntryPoint>]
let main _ =
    let blades = [
        0b10uy, 10f
        0b01uy, System.MathF.PI
        0uy, 1f
        0b11uy, 99f
    ]
    do testMultivectorOfBlades blades

    0