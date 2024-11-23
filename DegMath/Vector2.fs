namespace DegMath

module Vector2 =
    let add (x1, y1) (x2, y2) = x1+x2, y1+y2 :float32*float32
    let sub (x1, y1) (x2, y2) = x1-x2, y1-y2 :float32*float32
    let neg (x, y) = -x, -y :float32*float32
    
    let scale s (x,y) = s*x, s*y :float32*float32

    let dot (x1, y1) (x2, y2) = x1*x2 + y1*y2 :float32

    let magSqr (x,y) = (x*x + y*y) :float32
    let mag = magSqr >> System.MathF.Sqrt

    let normalize v =
        let vmag = mag v
        if (vmag = 0f) then (v, 0f) else
        scale (1f/vmag) v, vmag