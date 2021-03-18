iterator_field(angspe::FieldAngularSpectrumScalar) = angspe.e_SXY
iterator_index(angspe::FieldAngularSpectrumScalar) = 1:length(angspe.e_SXY)

iterator_field(space::FieldSpaceScalar) = space.e_SXY
iterator_index(space::FieldSpaceScalar) = 1:length(space.e_SXY)

iterator_index_field(space::AbstractFieldMonochromatic) = (iterator_index(space), iterator_field(space))
