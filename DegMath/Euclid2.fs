namespace DegMath
open Clifford

module Euclid2 =

    let private PGA2 = Cl(2,0,1)

    let private point (x, y) = [0b011uy, 1f; 0b110uy, x; 0b101uy, y]
    let private getPoint m_p =
        let rec loop (y0, x0, xy) = function
        | [] when (xy = 0f) -> 
            (0f,0f)
        | [] -> 
            let s = 1f/xy
            (y0*s, x0*s)

        | (0b110uy, mag) :: t -> loop (y0+mag, x0, xy) t
        | (0b101uy, mag) :: t -> loop (y0, x0+mag, xy) t
        | (0b011uy, mag) :: t -> loop (y0, x0, xy+mag) t

        | _ :: t -> loop (y0, x0, xy) t
        loop (0f, 0f, 0f) m_p

    type Transform = 
        private |Transform of Multivector

        member this.transform = 
            match this with
            Transform t -> point >> PGA2.Sandwich t >> getPoint

        member this.inverse =
            match this with 
            Transform t -> t |> PGA2.VersorInv |> Transform

    let identity = Transform [0uy, 1f]

    //exponentiating infinite points
    let translate (dx, dy) = 
        Transform [0uy, 1f; 0b101uy, 0.5f * dx; 0b110uy, -0.5f * dy]

    //exponentiating finite points
    let rotate invariantSubspace angle =
        let halfAngle = angle * 0.5f
        let axis = point invariantSubspace
        Transform ((0uy, cos halfAngle) :: PGA2.Scale (sin halfAngle) axis)

    let precedes b a =
        match a,b with
        Transform t1, Transform t2 -> 
            t1 
            |> PGA2.Mul t2 
            |> (PGA2.Normalize >> fst)
            |> Transform

    let succeeds b a = precedes a b