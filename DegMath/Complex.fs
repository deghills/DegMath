namespace DegMath

open System

type private Complex = float32 * float32

module Complex =
    let getReal (c: Complex) = fst c
    let getIm (c: Complex) = snd c

    let zero = 0f, 0f

    let cnj (r, i) = r, -i
    let magSqr (r, i) = r * r + i * i
    let mag = magSqr >> MathF.Sqrt

    let add (r1, i1) (r2, i2) = r1 + r2, i1 + i2
    let sub (r1, i1) (r2, i2) = r1 - r2, i1 - i2
    let mul (r1, i1) (r2, i2) = (r1 * r2) - (i1 * i2), (r1 * i2) + (i1 * r2)
    let scale s = mul (s, 0f)
    let inverseScale s = mul (1f/s, 0f)

    let inv c =
        let numerator = cnj c
        let denominator = 1f / magSqr c, 0f
        mul numerator denominator

    let div a b = mul a (inv b)

    let normalize c =
        let cmag = mag c

        if (cmag = 0f) 
        then failwith "DIVIDE BY ZERO EXCEPTION" 
        else inverseScale cmag c, cmag

    let sqrt c = 
        let chat, cmag = normalize c        
        let halfangle = chat |> add (1f, 0f) |> normalize |> fst
        halfangle |> (MathF.Sqrt >> scale) cmag

    let EulersForm x = MathF.Cos(x), MathF.Sin(x)